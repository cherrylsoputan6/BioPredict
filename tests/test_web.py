import unittest
import json
from http.client import HTTPConnection
from http.server import ThreadingHTTPServer
from threading import Event, Thread
from web_app import PredictionJobs, handler


class WebJobTests(unittest.TestCase):
    def test_http_routes_and_origin_checks(self):
        jobs = PredictionJobs(lambda g,v: dict(gene=g, variant=v, score_label_1=.2))
        server = ThreadingHTTPServer(('127.0.0.1', 0), handler(jobs))
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        connection = HTTPConnection('127.0.0.1', server.server_port, timeout=5)
        try:
            connection.request('GET','/')
            response = connection.getresponse()
            self.assertEqual(response.status,200)
            self.assertIn(b'Research', response.read())
            connection.request('GET','/data/conservation_model/model.joblib')
            response = connection.getresponse()
            self.assertEqual(response.status,404)
            response.read()
            body = json.dumps(dict(gene='EXAMPLE',variant='R2H'))
            connection.request('POST','/api/predict',body,{'Content-Type':'application/json','Origin':'https://example.com'})
            response = connection.getresponse()
            self.assertEqual(response.status,403)
            response.read()
            connection.request('POST','/api/predict',body,{'Content-Type':'application/json','Origin':f'http://127.0.0.1:{server.server_port}'})
            response = connection.getresponse()
            self.assertEqual(response.status,202)
            identity = json.loads(response.read())['job']
            jobs.jobs[identity].result(timeout=5)
            connection.request('GET','/api/jobs/'+identity)
            response = connection.getresponse()
            self.assertEqual(json.loads(response.read())['status'],'complete')
        finally:
            connection.close()
            server.shutdown()
            server.server_close()
            thread.join(timeout=5)
            jobs.close()

    def test_invalid_input_never_invokes_prediction(self):
        called = []
        jobs = PredictionJobs(lambda g,v: called.append((g,v)))
        try:
            for gene, variant in [('bad symbol','R2H'), ('EXAMPLE','M1I'), ('EXAMPLE','R2*')]:
                with self.assertRaises(ValueError):
                    jobs.submit(gene, variant)
            self.assertEqual(called, [])
        finally:
            jobs.close()

    def test_single_job_limit_and_failure_are_reported(self):
        waiting, release = Event(), Event()
        def predictor(gene, variant):
            waiting.set()
            release.wait(timeout=5)
            raise ValueError('Coverage missing')
        jobs = PredictionJobs(predictor)
        try:
            identity = jobs.submit('EXAMPLE','R2H')
            self.assertTrue(waiting.wait(timeout=5))
            self.assertEqual(jobs.status(identity)['status'],'running')
            with self.assertRaises(RuntimeError):
                jobs.submit('OTHER','R3H')
            release.set()
            try:
                jobs.jobs[identity].result(timeout=5)
            except ValueError:
                pass
            self.assertEqual(jobs.status(identity),dict(status='unavailable',error='Coverage missing'))
            self.assertIsNone(jobs.status('unknown'))
        finally:
            release.set()
            jobs.close()
