"""Local-only BioPredict research interface; frozen eight-feature inference."""
import json
import re
import secrets
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Lock
from predict_variant import parse_missense, predict

ROOT = Path(__file__).resolve().parent


class PredictionJobs:
    def __init__(self, predictor=None):
        self.pool = ThreadPoolExecutor(max_workers=1)
        self.lock = Lock()
        self.jobs = {}
        self.predictor = predictor or self._predict

    @staticmethod
    def _predict(gene, variant):
        return predict(gene, variant, ROOT / 'data/conservation_model/model.joblib',
                       ROOT / 'data/sequence_cache', ROOT / 'data/ortholog_cache', timeout=30)

    def submit(self, gene, variant):
        if not isinstance(gene, str) or not re.fullmatch(r'[A-Za-z0-9_-]{1,40}', gene):
            raise ValueError('Enter a gene symbol without spaces')
        if not isinstance(variant, str) or len(variant) > 20:
            raise ValueError('Enter a one-letter protein substitution')
        parse_missense(variant)
        with self.lock:
            if any(not job.done() for job in self.jobs.values()):
                raise RuntimeError('A prediction is already running; wait for it to finish')
            if len(self.jobs) >= 20:
                self.jobs.pop(next(iter(self.jobs)))
            identifier = secrets.token_urlsafe(18)
            self.jobs[identifier] = self.pool.submit(self.predictor, gene, variant)
            return identifier

    def status(self, identifier):
        with self.lock:
            job = self.jobs.get(identifier)
        if job is None:
            return None
        if not job.done():
            return dict(status='running')
        try:
            return dict(status='complete', result=job.result())
        except Exception as error:
            return dict(status='unavailable', error=str(error))

    def close(self):
        self.pool.shutdown(wait=True)


def handler(jobs):
    class Handler(BaseHTTPRequestHandler):
        def send(self, status, body, content_type='application/json'):
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(status)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.end_headers()
            self.wfile.write(data)

        def valid_host(self):
            port = self.server.server_port
            return self.headers.get('Host') in (f'127.0.0.1:{port}', f'localhost:{port}')

        def do_GET(self):
            if not self.valid_host():
                return self.send(403, dict(error='Local access only'))
            if self.path == '/':
                return self.send(200, (ROOT / 'web/index.html').read_bytes(), 'text/html; charset=utf-8')
            if self.path.startswith('/api/jobs/'):
                result = jobs.status(self.path.removeprefix('/api/jobs/'))
                return self.send(200 if result else 404, result or dict(error='Unknown job'))
            return self.send(404, dict(error='Not found'))

        def do_POST(self):
            origin = self.headers.get('Origin')
            allowed_origins = (f'http://127.0.0.1:{self.server.server_port}', f'http://localhost:{self.server.server_port}')
            if not self.valid_host() or origin not in allowed_origins:
                return self.send(403, dict(error='Local same-origin access required'))
            if self.path != '/api/predict':
                return self.send(404, dict(error='Not found'))
            try:
                length = int(self.headers.get('Content-Length', '0'))
                if not 0 < length <= 1024 or self.headers.get('Content-Type') != 'application/json':
                    raise ValueError('Expected a small JSON request')
                request = json.loads(self.rfile.read(length))
                if not isinstance(request, dict):
                    raise ValueError('Expected gene and variant fields')
                identifier = jobs.submit(request.get('gene'), request.get('variant'))
                return self.send(202, dict(job=identifier))
            except (ValueError, TypeError) as error:
                return self.send(400, dict(error=str(error)))
            except RuntimeError as error:
                return self.send(429, dict(error=str(error)))
    return Handler


def main():
    jobs = PredictionJobs()
    server = ThreadingHTTPServer(('127.0.0.1', 8000), handler(jobs))
    print('BioPredict research prototype: http://127.0.0.1:8000', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        jobs.close()


if __name__ == '__main__':
    main()
