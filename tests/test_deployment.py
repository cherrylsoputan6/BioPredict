import unittest
from unittest.mock import patch
from web_app import PredictionJobs
from deploy_app import create_app


class DeploymentTests(unittest.TestCase):
    def test_health_input_origin_and_expired_jobs(self):
        jobs = PredictionJobs(lambda g,v: dict(gene=g,variant=v,transcript='NM_123.1',
            model_class='benign',score_label_1=.2,threshold=.5,ortholog_species=['mouse'],
            interpretation='Uncalibrated',features={'relative_position':.5}))
        app = create_app(jobs)
        client = app.test_client()
        try:
            self.assertEqual(client.get('/health').status_code,200)
            with client.get('/') as home:
                self.assertEqual(home.status_code,200)
            self.assertEqual(client.get('/data/sequence_cache/a.gb').status_code,404)
            with patch.dict('os.environ',{'BIOPREDICT_PUBLIC_ORIGIN':'https://demo.example'}):
                self.assertEqual(client.post('/api/predict',json=dict(gene='TP53',variant='R175H'),headers={'Origin':'https://evil.example'}).status_code,403)
                self.assertEqual(client.post('/api/predict',json=dict(gene='TP53',variant='M1I'),headers={'Origin':'https://demo.example'}).status_code,400)
                response = client.post('/api/predict',json=dict(gene=' TP53 ',variant=' r175h '),headers={'Origin':'https://demo.example'})
                self.assertEqual(response.status_code,202)
                identity = response.json['job']
                jobs.jobs[identity].result(timeout=5)
                self.assertEqual(client.get('/api/jobs/'+identity).json['result']['variant'],'R175H')
                downloaded = client.get('/api/jobs/'+identity+'/report.json')
                self.assertEqual(downloaded.status_code,200)
                self.assertIn('attachment;',downloaded.headers['Content-Disposition'])
                self.assertEqual(downloaded.json['independent_validation'],'not_completed')
                readable = client.get('/api/jobs/'+identity+'/report.txt')
                self.assertIn(b'Uncalibrated label-1 score: 0.2',readable.data)
                self.assertEqual(client.get('/api/jobs/'+identity+'/report.exe').status_code,404)
                self.assertEqual(client.get('/api/jobs/missing/report.txt').status_code,409)
            self.assertEqual(client.get('/api/jobs/missing').status_code,404)
        finally:
            jobs.close()
