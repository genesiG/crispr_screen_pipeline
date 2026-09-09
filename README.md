# Pooled CRISPR Screen Analysis Pipeline

A config-driven, HPC-batch-scheduled pipeline for pooled CRISPR knockout screens: raw
single-end FASTQ through barcode/guide demultiplexing, sgRNA read counting, MAGeCK
statistical testing across arbitrary treatment/control comparisons, and CRISPRO-based
mapping of screen hits onto protein structure. It includes the same agentic orchestration
layer (Skills/, Specialists/, AGENTS.md) used across this pipeline family, adapted for a
CRISPR screen's step sequence.

This repository is a sanitized snapshot of a working research pipeline. Absolute paths, the
originating HPC username, and the real cell-line/condition labels used in one specific screen
have been replaced with placeholders, and directories holding raw sequencing data and real
sample/library metadata have been excluded.

## Architecture

```
Scripts/
├── config.py                  # Single source of truth: paths, barcode/guide column layout,
│                               # mismatch tolerances, and the TREATMENTS/CONTROLS comparison
│                               # matrix — every step reads from here instead of hardcoding values.
├── step_0_setup.py            # Directory scaffolding
├── step_0b_setup.sh           # Shell-side environment/setup companion
├── step_1_demultiplex.py      # Demultiplex pooled FASTQ by sample barcode and count sgRNA reads
├── step_1qc.py / step_1qc_report.py   # Demultiplexing QC + summary report
├── step_2_count.py            # Assemble the per-sample sgRNA count matrix
├── step_2qc.py                 # Guide-representation / coverage QC
├── step_3_mageck_test.py      # Iterates configured treatment/control pairs, writes and
│                               # submits one MAGeCK test batch job per comparison
├── step_4_run_crispro.py      # Submit CRISPRO batch jobs to map hit sgRNAs onto protein
│                               # structure/domain coordinates
├── step_5_color_pdb.py        # Render screen hits as per-residue coloring on a PDB structure
└── utils/guide_representation.R   # Guide-level representation/coverage statistics

Skills/, Specialists/, AGENTS.md   # Agentic orchestration layer (see below)
```

**Design principles**
- Paths are derived from `config.py`.
- Every step is per-sample and submits its own HPC batch job (`bsub`), writing its own
  `.batch`/`.log`/`.error` files so a failed sample can be re-run in isolation.
- `TREATMENTS`/`CONTROLS` are plain pattern lists matched against count-table column names, so
  adding a new comparison can be done by simply editing the `config.py` module
— `step_3_mageck_test.py` iterates them to submit one MAGeCK job per treatment/control pair.

## Agentic orchestration layer

`AGENTS.md`, `Specialists/`, and `Skills/` define the same multi-agent system (Orchestrator / QC
/ Engineer / Analyst) used across this pipeline family, here driving a CRISPR screen: submitting
demultiplexing and MAGeCK jobs, polling with exponential backoff, and handing off to QC after each
step — see the `cut_and_run_pipeline` and `rnaseq_pipeline` repos for the full description of the
Specialist roles.

## Usage

1. Create the conda environment: `conda env create -f environment.yml`
2. Edit `Scripts/config.py`: set `USERNAME`, `PROJECT_NAME`, the sgRNA library/barcode column
   layout, and your `TREATMENTS`/`CONTROLS` comparison lists.
3. Populate `Metadata/` with your gRNA library table(s) and barcode table (see
   `Scripts/guide_to_scripts.txt` for the expected layout).
4. Run steps in order, `step_0_setup.py` → `step_5_color_pdb.py`, or let the orchestrator agent
   drive the run — see `Skills/skill_orchestrator.md`.

## License

MIT — see [LICENSE](LICENSE).
