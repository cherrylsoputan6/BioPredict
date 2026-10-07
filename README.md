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

## Planned Features

- Automatic selection of canonical RefSeq transcripts
- FASTA file upload
- Protein analysis
- Interactive web application (Streamlit)
- AI-powered biological sequence prediction
- PDF report generation

## Author

Cherryl Soputan
