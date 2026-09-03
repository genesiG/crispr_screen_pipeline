#!/usr/bin/env python3
"""
step_1_demultiplex.py
======================
Barcode demultiplexing — splits the single multiplexed FASTQ file into one
compressed FASTQ per sample using the inline barcodes in Metadata/barcodes.csv.

Tool used : cutadapt  (multi-threaded, I/O-efficient, exact barcode matching)
LSF       : 1 job, 8 cores, 8 h wall time

Skip-if-exists:
    If every per-sample FASTQ already exists in Importable_Data/demultiplexed_fastqs/
    the script prints a message and exits without submitting a new job.

Run after this step completes:
    python Scripts/step_2_count.py

Inputs (from Scripts/config.py):
    config.FASTQ_PATH  — multiplexed .fastq.gz
    config.BARCODES    — Metadata/barcodes.csv  (columns configured via SAMPLE_NAME_COL and BARCODE_SEQUENCE_COL)

Outputs:
    Importable_Data/demultiplexed_fastqs/<sample_name>.fastq.gz  — one per sample
    Analysis_Data/counts/log/<PROJECT>_step1_demux.log           — LSF stdout
    Analysis_Data/counts/log/<PROJECT>_step1_demux.error         — LSF stderr
    Analysis_Data/counts/log/<PROJECT>_step1_cutadapt.log        — cutadapt output
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
        sys.exit(1)
    except FileNotFoundError:
        print(
            f"  [DRY-RUN] bsub not found — job script for {label}:\n{script}",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    print(f"Project : {config.PROJECT_NAME}")
    print(f"FASTQ   : {config.FASTQ_PATH}")
    print(f"Barcodes: {config.BARCODES}\n")

    # Ensure output directories exist
    demux_dir = os.path.join(config.IMPORTABLE_DATA, "demultiplexed_fastqs")
    log_dir  = os.path.join(demux_dir, "log")
    os.makedirs(log_dir,   exist_ok=True)
    os.makedirs(demux_dir, exist_ok=True)

    # ---- Read barcodes.csv --------------------------------------------------
    sample_fastqs: dict[str, str] = {}   # sample_name -> expected output path
    cutadapt_adapter_args: list[str] = []

    with open(config.BARCODES, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        if getattr(config, "HAS_HEADER", True):
            next(reader, None)
        for row in reader:
            if not row or len(row) < max(config.SAMPLE_NAME_COL, config.BARCODE_SEQUENCE_COL):
                continue
            name  = row[config.SAMPLE_NAME_COL - 1].strip()
            index = row[config.BARCODE_SEQUENCE_COL - 1].strip()
            out_fq = os.path.join(demux_dir, f"{name}.fastq.gz")
            sample_fastqs[name] = out_fq
            # -g NAME=SEQUENCE allows variable 5' stagger before the barcode
            # and trims both the preceding stagger and the barcode before writing.
            cutadapt_adapter_args.append(f"-g {name}={index}")

    print(f"Samples : {len(sample_fastqs)} barcodes loaded")
    for name, fq in sample_fastqs.items():
        status = "✓ exists" if os.path.exists(fq) else "○ pending (to be generated)"
        print(f"  {status}  {name}")

    # ---- Skip-if-exists guard -----------------------------------------------
    if all(os.path.exists(fq) for fq in sample_fastqs.values()):
        print(
            "\nAll per-sample FASTQs already exist — nothing to submit.\n"
            "Proceed directly with:  python Scripts/step_2_count.py"
        )
        return

    # ---- Build the bsub job script ------------------------------------------
    adapter_args_str      = " \\\n    ".join(cutadapt_adapter_args)
    demux_output_pattern  = os.path.join(demux_dir, "{name}.fastq.gz")

    log_out    = os.path.join(log_dir, f"{config.PROJECT_NAME}_step1_demux.log")
    log_err    = os.path.join(log_dir, f"{config.PROJECT_NAME}_step1_demux.error")
    cutadapt_log = os.path.join(log_dir, f"{config.PROJECT_NAME}_step1_cutadapt.log")

    search_in_header = getattr(config, "SEARCH_BARCODE_IN_HEADER", False)
    if search_in_header:
        mode_label = "header (@read_id...)"
        barcodes_dict_str = repr(sample_fastqs)
        demux_cmd = f"""python3 -c '
import gzip, sys
barcodes = {barcodes_dict_str}
handles = {{b: gzip.open(path, "wt", compresslevel=4) for b, path in barcodes.items()}}
counts = {{b: 0 for b in barcodes}}
with gzip.open("{config.FASTQ_PATH}", "rt") as f_in:
    while True:
        h = f_in.readline()
        if not h: break
        s = f_in.readline(); p = f_in.readline(); q = f_in.readline()
        for b, h_out in handles.items():
            if b in h:
                h_out.write(h + s + p + q)
                counts[b] += 1
                break
for h_out in handles.values(): h_out.close()
for b, c in counts.items(): print(f"  {{b}}: {{c}} reads")
' 2>&1 | tee "{cutadapt_log}" """
    else:
        mode_label = "read sequence (allowing stagger, exact matching)"
        demux_cmd = f"""cutadapt \\
    {adapter_args_str} \\
    -e 0 \\
    --discard-untrimmed \\
    --no-indels \\
    --cores 0 \\
    --output "{demux_output_pattern}" \\
    "{config.FASTQ_PATH}" \\
    2>&1 | tee "{cutadapt_log}" """

    job_script = f"""#!/bin/bash
#BSUB -J "step1_demux_{config.PROJECT_NAME}"
#BSUB -q normal
#BSUB -o "{log_out}"
#BSUB -e "{log_err}"
#BSUB -n 8
#BSUB -M 16000
#BSUB -R "rusage[mem=2000] span[hosts=1]"
#BSUB -W 8:00

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate crisprscreen

set -euo pipefail

echo "=== Step 1: barcode demultiplexing ({mode_label}) ==="
echo "Input  : {config.FASTQ_PATH}"
echo "Output : {demux_dir}"
echo "Started: $(date)"

{demux_cmd}

echo "Finished: $(date)"
echo "=== Per-sample read counts ==="
for fq in {demux_dir}/*.fastq.gz; do
    n=$(zcat "$fq" | awk 'NR%4==1' | wc -l)
    echo "  ${{fq##*/}}  $n reads"
done
echo ""
echo "Next step: python Scripts/step_2_count.py"
"""

    print("\nSubmitting barcode demultiplexing job (cutadapt)...")
    _bsub(job_script, label="step1_demux")
    print(
        "\nMonitor with:  bjobs -u $USER\n"
        "When the job finishes, run:  python Scripts/step_2_count.py"
    )


if __name__ == "__main__":
    main()
