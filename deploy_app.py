"""Single-worker WSGI research service with bounded jobs and same-origin requests."""
import json
import os
import time
from collections import deque
from threading import Lock
from flask import Flask, jsonify, request, send_file
from werkzeug.middleware.proxy_fix import ProxyFix
from web_app import PredictionJobs, ROOT
from portable_model import load_portable
from reporting import report_bytes
from flask import Response


def create_app(jobs=None):
    app = Flask(__name__)
    app.wsgi_app = ProxyFix(app.wsgi_app, x_proto=1)
    app.config['MAX_CONTENT_LENGTH'] = 1024
    jobs = jobs or PredictionJobs()
    started, lock = deque(), Lock()

    @app.get('/')
    def home():
        return send_file(ROOT / 'web/index.html', max_age=0)

    @app.get('/health')
    def health():
        try:
            model = load_portable(ROOT / 'artifacts/baseline.json')
            return jsonify(status='ready', model_id=model['model_id'], independent_validation='not_completed')
        except (ValueError, OSError, KeyError):
            return jsonify(status='unavailable'), 503

    @app.post('/api/predict')
    def submit():
        expected = os.getenv('BIOPREDICT_PUBLIC_ORIGIN') or os.getenv('RENDER_EXTERNAL_URL')
        if not expected or request.headers.get('Origin') != expected.rstrip('/'):
            return jsonify(error='Configured same-origin access required'), 403
        body = request.get_json(silent=True)
        if not isinstance(body, dict):
            return jsonify(error='Enter a gene symbol and protein substitution'), 400
        try:
            with lock:
                now = time.monotonic()
                while started and now - started[0] >= 3600:
                    started.popleft()
                if len(started) >= 20:
                    return jsonify(error='Research service hourly capacity reached. Try later.'), 429
                identifier = jobs.submit(body.get('gene'), body.get('variant'))
                started.append(now)
            return jsonify(job=identifier), 202
        except ValueError as error:
            return jsonify(error=str(error)), 400
        except RuntimeError as error:
            return jsonify(error=str(error)), 429

    @app.get('/api/jobs/<identifier>')
    def status(identifier):
        result = jobs.status(identifier)
        if result is None:
            return jsonify(error='Result expired or service restarted. Submit again.'), 404
        if result['status'] == 'unavailable':
            # Unexpected exceptions may expose local paths; keep the public error bounded.
            text = result['error']
            if any(word in text for word in ('timed out', 'getaddrinfo', 'HTTP Error')):
                result = dict(status='unavailable', error='Sequence provider unavailable. Retry later; no score was issued.')
            elif '/' in text or '\\' in text:
                result = dict(status='unavailable', error='Model or sequence data unavailable; no score was issued.')
        return jsonify(result)

    @app.get('/api/jobs/<identifier>/report.<kind>')
    def download(identifier, kind):
        if kind not in ('json', 'txt'):
            return jsonify(error='Unknown report format'), 404
        result = jobs.status(identifier)
        if result is None or result['status'] != 'complete':
            return jsonify(error='Completed report unavailable or expired'), 409
        data, content_type = report_bytes(result['result'], kind)
        response = Response(data, content_type=content_type)
        response.headers['Content-Disposition'] = f'attachment; filename="BioPredict-report.{kind}"'
        return response

    @app.after_request
    def headers(response):
        response.headers['Cache-Control'] = 'no-store'
        response.headers['X-Content-Type-Options'] = 'nosniff'
        response.headers['X-Frame-Options'] = 'DENY'
        response.headers['Referrer-Policy'] = 'no-referrer'
        return response

    return app


app = create_app()
