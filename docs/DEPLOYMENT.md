# Public research-demo deployment

## Current deployment attempt

On 2026-10-17, the authenticated Render new-service form showed the connected
BioPredict repository. The $0/month plan, Python 3.12.14, deployment requirements,
single-worker start command, `/health`, and auto-deploy Off were configured.
Submitting the form triggered an Add Card verification dialog stating a temporary
$1 USD authorization. No payment details were entered. The owner must handle
this billing verification or choose a different host before public deployment
can continue. No successful deployment or public URL has been verified.

The repository includes a [Render Blueprint](https://render.com/docs/blueprint-spec)
for a Python free-tier research demonstration. Hosting requires a Render account
authorized to connect this GitHub repository; no Render credentials are available
in the current workspace. Public deployment has **not** been performed or verified.
The browser check reached Render's sign-in page; an authenticated account session
is required. Local adapter and report endpoints have been tested (43 tests pass),
and the browser downloaded a real JSON attachment. A public HTTPS deployment and
Linux dependency build still require verification after account connection.

1. Sign in to Render and connect `cherrylsoputan6/BioPredict`.
2. Create a Blueprint from the tracked `render.yaml`; confirm the free plan and
   set `NCBI_EMAIL` to the contact used for NCBI requests. No paid resource is declared.
3. Build with `requirements-deploy.txt` and start the single Gunicorn worker.
   Render supplies `PORT` and `RENDER_EXTERNAL_URL` for health/origin checks.
4. Verify `/health` is ready, homepage rendering, invalid-input rejection, a
   successful sequence/coverage demonstration, refusal behavior and both report
   downloads on the actual HTTPS URL before recording deployment success.

The public model is a data-only JSON export verified against all 679 frozen
training predictions: maximum numerical difference 1.11e-16. It is the same
eight-feature baseline, not a newly trained model; no untrusted pickle is loaded
by the public service. Source artifact hash and frozen configuration are retained.

[Free Render services](https://render.com/docs/free) sleep on idle and lose local
filesystem changes on restart/redeploy. Caches are regenerable and reports download
in the browser. Jobs may disappear after restart; the UI reports this rather than
pretending a result completed. One worker, one active prediction, a 20-job/hour
global limit and a 20-result memory bound keep the prototype's upstream load small.
These limits reset on restart; this is not a hardened multi-user production service.
Public endpoints never serve arbitrary files, training rows or the binary artifact.

Users submit a gene symbol and protein substitution. NCBI receives the gene/
transcript request; Ensembl receives the gene ortholog request. No patient records
or identifiers should be submitted. Results are temporary server memory and
download locally; platform/upstream logs may still exist. The interface is a
research demonstration with no independent clinical validation, no disease-risk
probability claims and no medical decisions.

## Completion criteria

- Scientific: new cohort provenance/lock, leakage screens and refusal accounting;
  prespecified evaluation completed, or explicitly marked unavailable.
- Application: inputs, explanations, refusals, provenance, both downloads and
  accessibility/usability checked on the deployed URL.
- Operational: reproducible artifact/dependencies, public HTTPS health verification,
  bounded requests, documented ephemeral-state/privacy limitations.
- Project: complete test suite passes, source commits pushed, PROJECT_STATUS current.

The source is deployment-ready for testing; the account connection and actual
public verification remain blocked pending hosting access. Independent scientific
validation remains incomplete and must not be implied by deployment readiness.
