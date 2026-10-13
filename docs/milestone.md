# BioPredict research prediction milestone

BioPredict now connects a gene/variant input to MANE sequence validation,
eight numerical features, and a saved logistic-regression pipeline. The existing
interactive sequence-analysis app remains available through `python app.py`.

## Dataset and evaluation

- 1,000 selected candidates yielded 971 sequence-validated substitutions.
- Gene-disjoint split: 778 training and 193 test variants.
- Archive ortholog coverage: 679 training and 176 test variants available.
- Five-fold training-only CV on identical eligible rows: six features 63.3%
  accuracy versus eight features 82.6%; ROC AUC 0.689 versus 0.882.
- Frozen eight-feature model test accuracy: 144/176 (81.8%), ROC AUC 0.885.
- Test confusion matrix: 53 benign correct, 27 benign flagged as pathogenic,
  5 pathogenic missed, 91 pathogenic correct.

Coverage differs by label, and unavailable rows are excluded rather than imputed.
The test split was previously evaluated for the original three-feature baseline;
it is not an entirely new independent dataset. Model settings were selected
through training-only analysis. Scores are uncalibrated and are not clinical risk.

## Working command demonstration

The user ran `python predict_variant.py TP53 R175H` successfully. The command
validated NM_000546.6, obtained four qualifying species, and emitted model label 1
with score 0.7221624865 at the fixed 0.5 threshold. All four retained the reference
arginine; none had histidine at that aligned site. This demonstrates operation,
not independent performance: demonstration variants can overlap development data.

## Reproduce and inspect

Install `requirements.txt` in the same Python environment used to run commands.
The locally saved model requires its recorded scikit-learn version. Run
`python -m unittest discover -s tests -v` for offline behavior tests.
`python demo_checks.py` replays three benign development examples and reference/
position refusal checks using existing local caches, without network requests.
Its report is saved under `data/` and deliberately does not tune the model.

Keep source data, source checksums, split summaries, conservation caches, frozen
configuration, CV reports, evaluation metrics and the saved model locally.
Generated data and models are ignored by Git. The saved pipeline should only be
loaded from a trusted local source.

## Cached demonstration results

The three benign-labeled development examples replayed from cached protein and
ortholog records produced the following fixed-model outputs:

| Variant | Dataset label | Model label | Score for label 1 |
| --- | --- | --- | --- |
| ARMH4 Q542E | Benign/likely benign | Benign | 0.02916 |
| SLITRK1 S330A | Benign/likely benign | Benign | 0.000266 |
| BCL11B F623L | Benign | Pathogenic | 0.71151 |

BCL11B F623L is a false positive relative to the dataset label, even though it was
included in development data. These outputs demonstrate behavior rather than
estimate accuracy. Six cached checks deliberately using incorrect reference
residues or positions beyond the protein length were all refused. The predictor
also has a test verifying it never calls model scoring after a reference failure.
The complete offline suite now contains 31 passing tests. Detailed demonstration
results remain in the ignored `data/demo_checks.json` artifact.

## Next research work

Training-only out-of-fold diagnostics found 92 false positives among 320 benign
variants and 26 false negatives among 359 pathogenic variants. Brier score was
0.12993 versus 0.24918 for a fold-training prevalence predictor. A descriptive
feature-group audit found 70/71 benign variants with complete reference agreement
across observed species were flagged, versus 22/249 with partial/no agreement.
Additionally, 91/145 benign variants with no observed alternate residue were
flagged, versus 1/175 where the alternate occurred in at least one observed
species. These overlapping groups suggest the model relies heavily on ortholog
agreement. They do not establish a causal effect or justify test-based tuning.

Preserve this configuration as a baseline. Independent evaluation on new genes
and a larger, carefully defined cohort is needed before claiming generalization.
Further work should characterize false positives, prediction coverage and score
calibration using training-only validation rather than tuning on this test set.
