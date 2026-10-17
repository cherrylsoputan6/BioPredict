# Prespecified independent-evaluation protocol

Status: no qualifying independent clinical cohort has been collected or evaluated.
The existing 193-row split is historically evaluated; it is not an untouched test.
The frozen eight-feature model and threshold 0.5 remain the sole baseline.

The original-source registry has been built: 4,646,126 unique VariationIDs from
9,223,698 rows and 1,016 excluded GeneIDs. The tracked manifest records checksums
and the historical-cohort extension; the SQLite database stays local. No future
cohort has passed the gate. Source publication date must still be verified before
a temporal claim; a local filename or file modification date is not sufficient.

## Source investigation

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/release_cycle/) provides archived
monthly releases and nonarchived weekly updates. A later *archived* release with
verified URL, publication date and checksum can support prospective assessment.
A record's LastEvaluated date is not its first appearance. Merely downloading a
newer file, or using changed classifications for existing IDs, is not sufficient.

[ProteinGym](https://proteingym.org/) distinguishes functional-assay DMS benchmarks
from clinical-annotation benchmarks. Its processed ClinVar/gnomAD clinical labels
may overlap our underlying source; they cannot automatically be called external
validation. DMS scores can assess functional ranking separately, with assay-specific
definitions and eligible human proteins, but must not be relabeled as clinical
pathogenic/benign without justified mappings. No DMS or external clinical results
have been fabricated or claimed.

## Cohort lock and exclusion rules

1. Freeze the complete original ClinVar source, source checksum and **all** VariationIDs,
   not just the sampled training IDs. Keep an SQLite registry via
   `python cohort_registry.py data/variant_summary.txt.gz data/development_candidates.csv data/prospective_registry`.
   The development file must combine all historically inspected candidate cohorts
   and demonstration genes. The gate records unmatched IDs and must fail closed.
2. Select a strictly later archived release. Confirm it predates all cohort scoring
   and preserve the source URI/date/checksum. Exclude every VariationID in the
   baseline registry, regardless of later label changes.
3. Require a resolved stable NCBI GeneID absent from **all** development/evaluation
   cohorts, including historical pilot, test and demonstration work. Gene aliases
   must map to the same excluded identity. Novel genes alone do not prove independence.
4. Require human GRCh38 simple missense substitutions, position >1, the same four
   accepted classes/combined labels and at least multiple submitters without
   conflicts, expert-panel or practice-guideline review. Record the narrower
   evidence population; preserve exclusions and do not change labels after scoring.
5. Validate exact MANE protein/reference and fixed archive ortholog coverage. Compare
   full sequences against every development protein; exclude exact identity. Audit
   homologous families before locking a stronger family-independent cohort. A
   prespecified similarity method/threshold is required before outcome inspection;
   absence of exact identity does not imply absence of family overlap.
6. Lock candidate IDs, protein hashes, exclusions and cohort checksum **before**
   model predictions. Keep all eligible samples; do not rebalance or select easier
   variants based on scores. Require at least 50 evaluable variants of each class
   across at least 20 stable gene groups per class for a planned primary analysis;
   below this, report feasibility/descriptive results only. These are operational
   minimums, not a statistical power calculation.

## Prespecified outcomes

Primary: balanced accuracy at 0.5, with gene-cluster bootstrap uncertainty.
Secondary: ROC AUC, precision/recall/F1 per class, confusion matrix, Brier/log loss,
fixed five-bin reliability diagram, calibration intercept/slope with uncertainty,
and performance by review status and gene (with denominators and sparse-group
warnings). Report attempted, refused and evaluable counts by class, plus all
refusal reasons. Feature-envelope excursions are descriptive OOD flags, not a
validated OOD detector. No test-based threshold selection or calibration fitting.

## Current limitation

No later archived snapshot and verified first-appearance ledger have yet been
acquired. A future snapshot is not guaranteed to yield enough new, gene-disjoint,
protein-family-separated examples. External functional assays answer a different
question. Existing performance remains internal/exploratory; publication as a
research demonstration must display this limitation prominently.
