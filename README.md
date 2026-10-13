# 🧬 BioPredict

BioPredict is a Python-based bioinformatics toolkit that retrieves and analyzes biological sequences from the NCBI database.

## Features

- 🔍 Search human genes by name
- 🌐 Download DNA sequences from NCBI
- 🧬 Calculate GC content
- 🔢 Count nucleotides (A, T, G, C)
- 🔄 Generate reverse complements
- 🧫 Translate DNA into protein sequences

## Technologies

- Python
- Biopython
- NCBI Entrez API
- Git
- GitHub

## Current Status

🚧 In development

## ClinVar dataset candidates

The existing interactive pipeline runs with `python app.py`. Install its
dependency with `python -m pip install -r requirements.txt`.

The next phase starts with a separate, offline dataset builder:

1. Download the official human ClinVar
   [variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
   into `data/`. NCBI's
   [file guide](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/)
   describes the selected metadata and monthly release schedule.
2. Run `python dataset.py data/variant_summary.txt.gz data/clinvar_candidates.csv`.
3. Save the printed counts with your download date and source file checksum
   so you can reproduce the dataset. Data files are excluded from Git.

This reads the compressed file directly, keeps GRCh38 records with a single
named gene and a versioned RefSeq transcript, and recognizes simple three-letter
protein substitutions in the ClinVar Name field. It accepts Benign, Likely
benign, Pathogenic, Likely pathogenic, and the corresponding combined labels.
Uncertain, conflicting, and other classifications are excluded. Duplicate
Variation IDs are removed. Review status, identifiers, original name, and
evaluation date remain in the CSV. The command refuses to overwrite an existing
output. Exclusion counts show which records the conservative parser cannot use.

These are **labeled candidates, not a training-ready feature table**. The summary
does not contain protein sequences or every transcript consequence. Before
extracting features, validate each substitution against its exact transcript
version. For a MANE-only dataset, require an exact accession match with the MANE
transcript; do not apply another transcript's residue position to MANE. Then
remove duplicate protein substitutions, exclude any inconsistent labels across
those substitutions, and split by gene to assess generalization without gene
overlap. No model has been trained yet.

Run offline dataset tests with `python -m unittest discover -s tests -v`.

### Validate a first batch

Run `python validate_dataset.py data/clinvar_candidates_v2.csv data/validated_batch.csv --limit 100`.
This processes the first 100 candidates (an ordered debugging batch, not a
representative training sample). It requires internet and Biopython. Each gene's
current MANE accession is looked up once per run; downloaded human GenBank records
are cached in `data/sequence_cache/` for reuse. Version mismatches are excluded
rather than remapped. Complete CDS translations must match the protein annotation.
The output retains every processed candidate with its validation status and error
detail. Only `valid` rows receive the three variant features and protein length.
Network failures are recorded as `sequence_error`; retry using a new output file.
Keep the cache, original candidates, and download metadata. This batch still needs
duplicate/conflict handling and a gene-held-out split before training.

### Select a broader pilot

Run `python sample_dataset.py data/clinvar_candidates_v2.csv data/pilot_candidates.csv --per-label 50 --seed 42`.
The selector scans all candidates before sampling, excludes all position-one
substitutions, groups by gene + exact transcript version + substitution, drops
whole groups with both labels, and retains one ClinVar record per remaining group.
It chooses equal numbers from each class, taking one variant per gene per round.
The fixed seed reproduces the selection for the same input. A companion
`pilot_candidates.summary.json` records counts, settings, and the input checksum.
This deliberately balanced, gene-diverse pilot does not reflect natural prevalence
and is not a train/test split. Conflict checks cover only the candidate CSV and
exact transcript versions; they do not resolve all ClinVar disease interpretations.

Validate it with `python validate_dataset.py data/pilot_candidates.csv data/pilot_validated.csv --limit 100 --timeout 60`.
Both classes must be counted again after MANE and sequence filtering. Keep only
valid rows for subsequent preparation, and use gene-disjoint splits before modeling.

### Prepare the pilot split

Run `python prepare_dataset.py data/pilot_validated_v2.csv data/pilot_split`.
This writes `clean.csv`, `train.csv`, `test.csv`, and `summary.json` into a new
folder. It retains only valid rows, rejects invalid numerical features and
duplicate substitutions, and places about 20% of genes in the test set. Both
sets must contain both labels, with no shared genes. The summary records label
counts, source checksum, seed, and the three allowed model features. Identifiers,
gene names, and clinical metadata remain for auditing; do not use them as model
inputs. Fit any scaling and model only on training data. The small, deliberately
balanced pilot is for checking the workflow, not estimating clinical performance.

### First baseline model

Install updated dependencies with `python -m pip install -r requirements.txt`.
Run `python train_model.py data/pilot_split data/pilot_model`.
The fixed logistic regression uses only relative position, hydrophobicity change,
and weight change. A StandardScaler pipeline fits only the training rows, following
[scikit-learn's leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html).
The command rejects overlapping genes, invalid features, duplicates within each
split, missing classes, and existing output folders. It saves `model.joblib`,
`metrics.json`, and `test_predictions.csv`. Accuracy, balanced accuracy, precision,
recall, F1, ROC AUC, and the confusion matrix are compared with a majority-class
baseline. The decision threshold is fixed at 0.5; there is no test-based tuning.
Scores from this deliberately balanced pilot are not calibrated clinical risks.
Load only trusted local joblib models; serialized models can execute code.

### Recover network failures

Run `python retry_validation.py extract data/larger_validated.csv data/larger_retry_candidates.csv`
to extract only recognized temporary connection errors. Validate that candidate
file using the normal validator with a limit large enough to include all rows.
After it finishes, run `python retry_validation.py merge data/larger_validated.csv data/larger_retried.csv data/larger_validated_merged.csv`.
The merge requires every retry candidate, rejects changed provenance, preserves
successful original rows and annotation exclusions, and refuses to overwrite files.
If connection errors remain, repeat extraction from the merged result with new
filenames. DNS failures (`getaddrinfo failed`) indicate hostname resolution failed;
retry after connectivity recovers rather than treating these as biological failures.

### Training-only cross-validation

Run `python cross_validate_model.py data/larger_split/train.csv data/larger_cv.json`.
Five stratified group folds hold out whole genes. Each fold fits its scaler and
fixed logistic regression only on that fold's training rows. The report compares
the baseline with a majority-class predictor, gives fold means and standard
deviations, and reports predictions made while each row was held out. Standard
deviations describe fold variability, not confidence intervals. Keep test.csv
out of cross-validation and feature selection. This command does not tune or
overwrite the saved model.

### Extended substitution descriptors

Run `python cross_validate_model.py data/larger_split/train.csv data/larger_cv_extended.json --extended`.
This uses the same folds and fixed model settings while adding absolute
hydrophobicity change, absolute weight change, and BLOSUM62 substitution score.
Absolute changes describe magnitude regardless of direction. The score comes from
[Biopython's bundled BLOSUM62 matrix](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html#loading-predefined-substitution-matrices).
It is a general substitution descriptor, not a site-specific conservation score
or a pretrained pathogenicity predictor. Feature computation checks the stored
deltas against the substitution and uses no labels or learned statistics. The
existing app and three-feature training command retain their original behavior.
Compare against the original CV report; do not tune using test-set performance.

### Audit label evidence and feature groups

Run `python audit_features.py data/larger_split/train.csv data/feature_audit`.
The audit counts clinical classifications and review statuses by label, then
compares original features, absolute changes added, BLOSUM62 added, and all six
features using identical training-only group folds and fixed settings. Each
comparison and a combined summary are saved in the new folder. Review status
describes submitted evidence strength; this metadata audit cannot establish
whether a label is correct. Repeated feature experiments can overfit these CV
folds, so use a preplanned selection and new independent data for final assessment.

### Review-status sensitivity check

Run `python evidence_sensitivity.py data/larger_split/train.csv data/evidence_sensitivity`.
This makes a separate training CSV excluding the exact review status
`no assertion criteria provided` and evaluates the original and extended feature
sets with gene-held-out CV. It records class counts before and after filtering and
compares both models with a majority baseline. Filtering changes the label and
gene composition and fold assignments, so a score difference is not proof of
label correctness or of improved generalization. No test rows or saved models
are modified. Focus on balanced accuracy and ROC AUC when comparing populations
with different class proportions.

## Planned Features

- Automatic selection of canonical RefSeq transcripts
- FASTA file upload
- Protein analysis
- Interactive web application (Streamlit)
- AI-powered biological sequence prediction
- PDF report generation

## Ortholog agreement pilot

Run `python conservation.py data/larger_split/train.csv data/conservation_pilot.csv --limit 10 --timeout 30`.
The script uses [Ensembl ortholog alignments](https://rest.ensembl.org/documentation/info/homology_symbol)
and the cached validated MANE proteins. It maps each human residue through gaps,
requires an exact full-length human protein match, and includes only one-to-one
orthologs from mouse, dog, cow, chicken, and zebrafish. At least three species must
have canonical residues at that column before a score is emitted. The resulting
fractions describe reference/alternate agreement, not phylogenetically weighted
conservation or clinical risk. Missing coverage is recorded, never replaced by
zero. Cached responses preserve URLs, download times and response hashes.
This initial coverage check uses only training variants; no model incorporates
these descriptors yet. Any later comparison must use identical eligible rows
and folds, because missing ortholog coverage can bias the sample.
For uncached genes, the fetcher requests the five species separately and saves
each completed response, reducing response size and preserving partial progress.
Existing full-response caches are still supported. A species request failure
marks the row as an error rather than silently computing a partial-panel score.
`--server oct2024` uses the official October 2024 Ensembl archive, with responses
in a separate cache subfolder. Archived annotations may differ; the exact MANE
protein check remains mandatory. Do not silently mix releases in a modeling
dataset. A backup-server probe must establish coverage before choosing a release
for the full training cohort.

Compare with `python compare_conservation.py data/larger_split/train.csv data/conservation_archive_training.csv data/conservation_comparison`.
This verifies candidate identities and all original training values, selects
available rows in original training order, and evaluates six features versus six
plus reference/alternate ortholog agreement on identical gene-held-out folds.
Coverage and label counts are reported. Agreement descriptors are not weighted
for evolutionary relatedness; any observed benefit applies to the covered
training subset and needs independent evaluation. The saved model and test set
remain unchanged.

The eight-feature evaluation configuration is recorded in `model_config.json`.
Before evaluating it, collect the existing test split with
`python conservation.py data/larger_split/test.csv data/conservation_archive_test.csv --limit 193 --timeout 60 --server oct2024`.
Collection does not fit a model. Preserve the model settings and report missing
coverage by class. This test split has previously been used for the original
three-feature baseline; it is not a new independent evaluation dataset.

Evaluate with `python evaluate_conservation.py data/conservation_comparison/matched_training.csv data/larger_split/test.csv data/conservation_archive_test.csv data/conservation_model`.
The command checks the frozen training checksum and original test provenance,
rejects overlapping genes, fits the fixed six/eight-feature models on identical
training rows, and evaluates both on identical covered test rows. It saves the
eight-feature pipeline, prediction table, coverage breakdown and metrics. The
configuration is not tuned using test results. Performance excludes unavailable
test rows and does not establish calibrated clinical risk.

## Research prediction command

Run `python predict_variant.py TP53 R175H` with a gene and one-letter substitution.
The example is not hard-coded; any qualifying human MANE gene can be supplied.
The separate command preserves `app.py`, validates the reference amino acid,
requires exact archived alignment coverage from at least three species, and uses
the saved eight-feature pipeline with its training scaler and fixed threshold.
It refuses to predict when coverage is missing and checks the training library
version before using the model. Network failures remain possible; caches are reused.
Load only the trusted local model created by this project. Output is a research
classification and uncalibrated score, not disease risk or a clinical decision.
A demonstration variant may have appeared in development data and therefore
does not provide new evidence of generalization.

## Author

Cherryl Soputan
