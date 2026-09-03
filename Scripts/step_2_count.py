#!/usr/bin/env python3
"""
step_2_count.py
================
sgRNA counting — for every gRNA library pool CSV, submits one LSF job that
runs `mageck count` over the per-sample FASTQs produced by step_1_demultiplex.py.

Tool used : mageck count  (multi-threaded; one job per library pool)
LSF       : 1 job per pool, 8 cores, 4 h wall time

Prerequisites:
    Run step_1_demultiplex.py and wait for it to finish before this step.

Skip-if-exists:
    Libraries whose .count.txt output file already exists are skipped silently.

Inputs (from Scripts/config.py):
    config.LIBRARY         — list of Metadata/grna_library_pool*.csv
                             (columns configured via GENE_SYMBOL_COL, ID_COL, TARGET_SEQUENCE_COL)
    config.TARGET_MISMATCHES — mismatch tolerance (1 recommended; capped at 1)
    Importable_Data/demultiplexed_fastqs/<sample>.fastq.gz  — from step 1

Outputs (per library pool):
    Analysis_Data/counts/<PROJECT>_<lib>_<N>mm.count.txt              — raw counts
    Analysis_Data/counts/<PROJECT>_<lib>_<N>mm.count_normalized.txt   — normalised
    Analysis_Data/counts/<PROJECT>_<lib>_<N>mm.log                    — mageck log
    Analysis_Data/counts/log/<PROJECT>_<lib>_step2.log                — LSF stdout
    Analysis_Data/counts/log/<PROJECT>_<lib>_step2.error              — LSF stderr
    Analysis_Data/counts/log/<PROJECT>_<lib>_mageck_count.log         — mageck output
    Importable_Data/mageck_libraries/<lib>_mageck_library.txt         — converted library
"""

import csv
import os
import subprocess
import sys

# ---------------------------------------------------------------------------
# Import pipeline config (Scripts/ directory)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
import config


# ---------------------------------------------------------------------------
# Helper: submit a job script string via bsub
# ---------------------------------------------------------------------------
def _bsub(script: str, label: str) -> None:
    """Submit *script* to LSF via bsub, printing the assigned job ID."""
    try:
        result = subprocess.run(
            ["bsub"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        print(f"  [submitted] {label}: {result.stdout.strip()}")
    except subprocess.CalledProcessError as exc:
        print(f"  [ERROR] submission failed for {label}:\n{exc.stderr}", file=sys.stderr)
    except FileNotFoundError:
        print(
            f"  [DRY-RUN] bsub not found — job script for {label}:\n{script}",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Helper: convert pool CSV to MAGeCK library format
# ---------------------------------------------------------------------------
def _build_mageck_library(lib_csv: str, out_dir: str, lib_name: str) -> str:
    """
    Convert <Gene, gRNA ID, Target Sequence, Library> CSV to a MAGeCK-ready
    tab-separated file with columns <sgRNA_id  sequence  gene>.

    Returns the path to the converted file. Skip-if-exists: if the file
    already exists it is returned immediately without re-reading the CSV.
    """
    out_path = os.path.join(out_dir, f"{lib_name}_mageck_library.txt")
    if os.path.exists(out_path):
        return out_path

    with (
        open(lib_csv, newline="", encoding="utf-8-sig") as fh_in,
        open(out_path, "w", encoding="utf-8") as fh_out,
    ):
        reader = csv.reader(fh_in)
        # Skip header line if present
        if getattr(config, "HAS_HEADER", True):
            next(reader, None)
        fh_out.write("sgRNA_id\tsequence\tgene\n")
        for i, row in enumerate(reader, start=1):
            if not row or len(row) < max([c for c in [config.GENE_SYMBOL_COL, config.TARGET_SEQUENCE_COL] if c is not None]):
                continue
            # Strip Windows carriage returns from every field
            gene = row[config.GENE_SYMBOL_COL - 1].strip().strip("\r")
            if config.ID_COL is not None and len(row) >= config.ID_COL:
                gid = row[config.ID_COL - 1].strip().strip("\r")
            else:
                gid = f"{gene}_{i}"
            seq  = row[config.TARGET_SEQUENCE_COL - 1].strip().strip("\r")
            if not seq:
                continue
            fh_out.write(f"{gid}\t{seq}\t{gene}\n")

    return out_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"Project  : {config.PROJECT_NAME}")
    print(f"Libraries: {len(config.LIBRARY)} pool(s) found")

    if not config.LIBRARY:
        print("ERROR: no library files found in config.LIBRARY — check config.py.", file=sys.stderr)
        sys.exit(1)

    # ---- Mismatch setting ---------------------------------------------------
    mismatches = config.TARGET_MISMATCHES
    if mismatches > 1:
        print(
            f"  [WARNING] config.TARGET_MISMATCHES={mismatches}; "
            "capping at 1 (2 is non-standard for 20 bp guides)."
        )
        mismatches = 1
    mismatch_label = f"{mismatches}mm"
    print(f"Mismatches: {mismatches}\n")

    # ---- Locate per-sample FASTQs from step 1 -------------------------------
    demux_dir = os.path.join(config.IMPORTABLE_DATA, "demultiplexed_fastqs")
    sample_fastqs: dict[str, str] = {}

    with open(config.BARCODES, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        if getattr(config, "HAS_HEADER", True):
            next(reader, None)
        for row in reader:
            if not row or len(row) < config.SAMPLE_NAME_COL:
                continue
            name  = row[config.SAMPLE_NAME_COL - 1].strip()
            out_fq = os.path.join(demux_dir, f"{name}.fastq.gz")
            sample_fastqs[name] = out_fq

    # Check that step 1 has actually been run
    missing = [fq for fq in sample_fastqs.values() if not os.path.exists(fq)]
    if missing:
        print(
            "ERROR: the following per-sample FASTQs from step 1 do not exist:\n"
            + "\n".join(f"  {f}" for f in missing)
            + "\n\nRun step_1_demultiplex.py and wait for it to finish first.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Per-sample FASTQs : {len(sample_fastqs)} files found in {demux_dir}")

    # ---- Directories --------------------------------------------------------
    log_dir    = os.path.join(config.COUNTDIR, "log")
    lib_scratch = os.path.join(config.IMPORTABLE_DATA, "mageck_libraries")
    os.makedirs(log_dir,     exist_ok=True)
    os.makedirs(lib_scratch, exist_ok=True)

    # ---- Build sample-label and fastq argument strings (shared across pools) -
    # All samples are passed to every mageck count run; pools that used a
    # different library will simply have near-zero counts in their columns.
    sample_label_str = ",".join(sample_fastqs.keys())
    fastq_str        = " ".join(f'"{fq}"' for fq in sample_fastqs.values())

    # ---- One job per library pool -------------------------------------------
    submitted = 0
    skipped   = 0
    n_cores   = 8

    for lib_path in sorted(config.LIBRARY):
        lib_name = os.path.splitext(os.path.basename(lib_path))[0]

        output_prefix = os.path.join(
            config.COUNTDIR,
            f"{config.PROJECT_NAME}_{lib_name}_{mismatch_label}",
        )
        count_table = output_prefix + ".count.txt"

        # ---- Skip-if-exists guard -------------------------------------------
        if os.path.exists(count_table):
            print(f"  [skip] {lib_name}: {count_table} already exists.")
            skipped += 1
            continue

        # Convert pool CSV → MAGeCK library format (skip if already done)
        mageck_lib = _build_mageck_library(lib_path, lib_scratch, lib_name)

        log_out    = os.path.join(log_dir, f"{config.PROJECT_NAME}_{lib_name}_step2.log")
        log_err    = os.path.join(log_dir, f"{config.PROJECT_NAME}_{lib_name}_step2.error")
        mageck_log = os.path.join(log_dir, f"{config.PROJECT_NAME}_{lib_name}_mageck_count.log")

        job_script = f"""#!/bin/bash
#BSUB -J "step2_count_{lib_name}"
#BSUB -q normal
#BSUB -o "{log_out}"
#BSUB -e "{log_err}"
#BSUB -n {n_cores}
#BSUB -M 32000
#BSUB -R "rusage[mem=4000] span[hosts=1]"
#BSUB -W 4:00

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate crisprscreen

set -euo pipefail

echo "=== Step 2: mageck count for {lib_name} ==="
echo "Library: {mageck_lib}"
echo "Started: $(date)"

mageck count \\
    --list-seq      "{mageck_lib}" \\
    --fastq         {fastq_str} \\
    --sample-label  "{sample_label_str}" \\
    --output-prefix "{output_prefix}" \\
    --sgrna-len     20 \\
    --trim-5        35 \\
    --norm-method   median \\
    2>&1 | tee "{mageck_log}"

echo "Finished: $(date)"
echo "=== Output files ==="
ls -lh "{output_prefix}"*
"""

        print(f"  Submitting mageck count job for {lib_name}...")
        _bsub(job_script, label=f"step2_count_{lib_name}")
        submitted += 1

    # ---- Summary ------------------------------------------------------------
    print(f"\nDone.  Submitted: {submitted}  Skipped: {skipped}")
    if submitted:
        print("Monitor with:  bjobs -u $USER")


if __name__ == "__main__":
    main()
