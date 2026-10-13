# BioPredict project status

## Current decision

Keep the frozen eight-feature logistic-regression model as the research baseline.
The three local sequence-context features remain experimental. Do not select a
replacement based on the existing test split or conditional bootstrap intervals.
No final clinical performance or probability-calibration claim is justified.

## Scientific validation milestone

The corrected training-only audit uses 679 variants, 671 gene symbols and 670
biological groups. Group components connect shared NCBI GeneIDs and identical
protein sequences, as well as gene symbols. This discovered IFT54/TRAF3IP1 alias
leakage in the earlier symbol-only folds: both use GeneID 26146 and NM_015650.4.
The original results remain historical exploratory estimates. They are not
silently overwritten. Corrected runs have zero gene/identical-protein overlap.

Three fixed seeds (42, 17, 123) each produce five paired grouped folds. Scaling is
fitted separately inside each training fold. Baseline and context use identical
held-out rows and fixed C=1, lbfgs, 1,000 iterations, threshold 0.5. Labels,
review statuses, identifiers and gene names do not enter model inputs. Local
context is calculated from the validated reference protein, excluding the variant
site. Ortholog scores are label-independent, fixed-species agreement descriptors;
they are not independent phylogenetic observations or site-specific clinical scores.

For seed 42, context minus baseline differences with conditional paired
gene-cluster 95% percentile intervals are:

| Metric | Difference | Interval |
| --- | ---: | --- |
| Accuracy | +0.88 percentage points | −0.73 to +2.50 points |
| Balanced accuracy | +0.92 percentage points | −0.70 to +2.53 points |
| ROC AUC | +0.0151 | +0.0034 to +0.0273 |
| Brier loss | −0.00683 | −0.01253 to −0.00120 |

The 2,000-draw paired bootstrap conditions on already generated OOF predictions.
It does not refit models, capture all training uncertainty or account for prior
adaptive feature experiments. Repeated folds are dependent, not 15 independent
experiments. ROC AUC/Brier signals support further research, not replacement.
The shuffled-label diagnostic has AUC 0.519 and balanced accuracy 0.505. This is
consistent with near-chance performance but does not prove leakage is absent.

## Labels, dataset bias and conserved benign variants

Selection deliberately balances labels and spreads candidates across genes.
It does not reproduce clinical prevalence. Labels aggregate ClinVar submissions
across conditions, and transcript-specific duplicate/conflict checks do not resolve
all disease-context disagreements. A benign label is not a guarantee about every
phenotype. Review-status counts describe evidence, not verified correctness.

Of 71 benign examples at fully conserved observed positions, 38 have a single
submitter with criteria, 28 multiple submitters without conflicts, and five no
assertion criteria. Thus the earlier 70/71 false-positive pattern cannot simply
be attributed to missing assertion criteria. Five-species agreement can miss
functional, structural and disease-context distinctions. Gene-diverse sampling
also produces genes with mostly one class and may create gene/function confounding.
These are hypotheses; no label has been changed to improve accuracy.

The corrected audit provides per-review, per-class and per-gene error counts.
663 gene symbols have one variant and eight have two, so reliable per-gene AUCs
cannot be claimed. Expert-panel coverage is only one training row.

## Calibration and applicability

Reliability bins and Brier loss use gene-held-out training predictions. Brier
combines calibration and discrimination; it does not establish calibration alone.
No calibrator or threshold adjustment is deployed. Scores must not be interpreted
as disease probabilities in a population with different prevalence.

Unsupported/OOD domains include nonhuman variants, nonmissense or position-one
changes, non-MANE isoforms, noncanonical residues, insufficient ortholog coverage,
novel protein families, temporal database changes and clinical populations with
different prevalence/review composition. Missing conservation is more common
among benign variants. Current results apply to the covered selected cohort.
Gene separation does not establish family-level independence. The predictor
refuses unsupported sequence/coverage cases but has no validated general OOD
detector. No disease-risk or calibrated-confidence display should be added to UI.

## Test-set lifecycle

The existing 193-row split was already evaluated by the three-feature and
eight-feature models. It cannot be described as untouched or reused for model
selection. This milestone loads training rows only. Further test evaluation is
frozen. A new cohort disjoint from all development genes, checked by stable gene
identity and preferably protein family, must be reserved before final selection.
Its design and provenance must be fixed before performance is inspected.

## Reproducibility and files

Run `python scientific_validation.py data/conservation_comparison/matched_training.csv data/scientific_validation_new`.
Run `python -m unittest discover -s tests -v` after changes.
The tracked scientific summary is in `reports/scientific_validation_summary.json`;
full fold gene lists, OOF predictions and per-gene tables remain under ignored
`data/scientific_validation_corrected/`. Source checksums and seeds are recorded.
The earlier milestone and demonstration details remain in `docs/milestone.md`.
Source data, caches and trained models are retained locally, not pushed to GitHub.

## Web application progress

A local research-only interface is available through `python web_app.py` on
127.0.0.1:8000. It wraps the existing eight-feature predictor, validates inputs,
handles requests asynchronously with a single-job limit, and displays coverage/
refusal messages and an uncalibrated model score. It does not load the context
model, retrain, tune, upload user results or publish a site. Same-origin/Host
checks and strict routes prevent serving arbitrary local files. Public deployment
and clinical-risk claims remain out of scope until independent validation.

## Next development

1. Reserve a genuinely new validation cohort before further selection; audit gene
   aliases and protein-family similarity across its boundary.
2. Validate generalization and subgroup coverage with prespecified metrics.
3. Compare calibration only with nested gene-group validation or a separate
   calibration cohort; retain the baseline unless evidence supports a change.
4. Build a local research web interface around the existing predictor with explicit
   coverage/refusal messages and provenance, without clinical-risk claims.
5. Deploy only after model/data lifecycle and independent validation are reviewable.
