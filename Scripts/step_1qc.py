#!/usr/bin/env python3
"""
step_1qc.py
===========
Step 1 Demultiplexing QC — submit all independent bsub jobs and extract
demultiplexing metrics (no inter-step dependencies within this script).

Run this FIRST after step_1_demultiplex.py, wait for all FastQC batch jobs
to finish, and then run:
    python3 Scripts/step_1qc_report.py

Steps (each outputting to its own subdirectory inside config.QC_DIR):
  fastqc          → FastQC on every demultiplexed FASTQ file (bsub) -> qc/fastqc
  cutadapt_stats  → Extract read fate & uniformity from cutadapt logs -> qc/cutadapt_stats

Usage:
    python3 Scripts/step_1qc.py                   # run all steps
    python3 Scripts/step_1qc.py --only fastqc     # submit FastQC jobs only
    python3 Scripts/step_1qc.py --only cutadapt_stats # extract metrics only
"""

import argparse
import csv
import glob
import json
import math
import os
import re
import sys
import time

# ---------------------------------------------------------------------------
# Import pipeline config (Scripts/ directory)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
import config

# ---------------------------------------------------------------------------
# Step registry & Directory paths
# ---------------------------------------------------------------------------
QC_STEPS = [
    ("fastqc",         "FastQC on per-sample demultiplexed FASTQs (bsub) -> qc/fastqc"),
    ("cutadapt_stats", "Extract demultiplexing metrics from cutadapt log -> qc/cutadapt_stats"),
]

DEMUX_DIR   = os.path.join(config.IMPORTABLE_DATA, "demultiplexed_fastqs")
FASTQC_DIR  = os.path.join(config.QC_DIR, "fastqc")
METRICS_DIR = os.path.join(config.QC_DIR, "cutadapt_stats")


# ===========================================================================
# Helper: calculate Gini coefficient & uniformity metrics
# ===========================================================================
def calc_uniformity_metrics(counts: list[int]) -> dict:
    valid_counts = sorted([c for c in counts if c > 0])
    n = len(valid_counts)
    if n == 0 or sum(valid_counts) == 0:
        return {"gini": 0.0, "cv": 0.0, "min_max_fold": 0.0, "mean": 0.0, "median": 0.0}

    total = sum(valid_counts)
    mean  = total / n
    stdev = math.sqrt(sum((c - mean) ** 2 for c in valid_counts) / n) if n > 1 else 0.0
    cv    = stdev / mean if mean > 0 else 0.0

    # Gini index
    gini = sum((2 * i - n - 1) * c for i, c in enumerate(valid_counts, 1)) / (n * total)

    min_val = valid_counts[0]
    max_val = valid_counts[-1]
    median  = valid_counts[n // 2] if n % 2 != 0 else (valid_counts[n // 2 - 1] + valid_counts[n // 2]) / 2.0
    min_max_fold = (max_val / min_val) if min_val > 0 else float("inf")

    return {
        "gini": round(gini, 4),
        "cv": round(cv, 4),
        "min_max_fold": round(min_max_fold, 2),
        "mean": round(mean, 1),
        "median": round(median, 1),
        "min": min_val,
        "max": max_val,
        "total": total,
    }


# ===========================================================================
# fastqc — FastQC on per-sample demultiplexed FASTQs -> qc/fastqc
# ===========================================================================
def _submit_fastqc_job(sample_name: str, in_file: str, out_dir: str, log_dir: str) -> None:
    job_name = f"fastqc_{sample_name}_{config.PROJECT_NAME}"
    log_out  = os.path.join(log_dir, f"{sample_name}_fastqc.log")
    log_err  = os.path.join(log_dir, f"{sample_name}_fastqc.error")
    batch_f  = os.path.join(log_dir, f"{sample_name}_fastqc.batch")

    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J "{job_name}"
#BSUB -o "{log_out}"
#BSUB -e "{log_err}"
#BSUB -q normal
#BSUB -n 2
#BSUB -M 8000
#BSUB -R "rusage[mem=2000] span[hosts=1]"
#BSUB -W 1:00

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate crisprscreen

set -euo pipefail

echo "=== FastQC: {sample_name} ==="
echo "Input  : {in_file}"
echo "Output : {out_dir}"
echo "Started: $(date)"

fastqc -o "{out_dir}" --extract --threads 2 "{in_file}"

rm -f "{out_dir}"/*fastqc.zip

echo "Finished: $(date)"
"""
    with open(batch_f, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"  Submitting FastQC job for {sample_name}...")
    os.system(f"bsub < {batch_f}")
    time.sleep(0.2)


def run_fastqc() -> None:
    """Submit one FastQC bsub job per demultiplexed sample writing into FASTQC_DIR."""
    log_dir = os.path.join(FASTQC_DIR, "log")
    os.makedirs(FASTQC_DIR, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    samples: list[str] = []
    if not os.path.exists(config.BARCODES):
        print(f"  [ERROR] Barcodes file not found: {config.BARCODES}", file=sys.stderr)
        return

    with open(config.BARCODES, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        if getattr(config, "HAS_HEADER", True):
            next(reader, None)
        for row in reader:
            if not row or len(row) < config.SAMPLE_NAME_COL:
                continue
            name = row[config.SAMPLE_NAME_COL - 1].strip()
            if name:
                samples.append(name)

    print(f"  Found {len(samples)} sample(s) defined in {config.BARCODES}")
    submitted = 0
    skipped   = 0

    for sample_name in sorted(samples):
        in_file     = os.path.join(DEMUX_DIR, f"{sample_name}.fastq.gz")
        out_html    = os.path.join(FASTQC_DIR, f"{sample_name}_fastqc.html")
        out_dir_ext = os.path.join(FASTQC_DIR, f"{sample_name}_fastqc")
        out_zip     = os.path.join(FASTQC_DIR, f"{sample_name}_fastqc.zip")

        if os.path.exists(out_html) and (os.path.isdir(out_dir_ext) or os.path.exists(out_zip)):
            print(f"  [skip] {sample_name}: FastQC report and data already exist in {FASTQC_DIR}.")
            skipped += 1
            continue

        if not os.path.exists(in_file):
            print(f"  [skip] {sample_name}: input FASTQ ({in_file}) not found. Run step_1_demultiplex.py first.")
            skipped += 1
            continue

        _submit_fastqc_job(sample_name, in_file, FASTQC_DIR, log_dir)
        submitted += 1

    print(f"  [fastqc] Submitted: {submitted} | Skipped: {skipped} | Output Dir: {FASTQC_DIR}")


# ===========================================================================
# cutadapt_stats — extract detailed demultiplexing metrics -> qc/cutadapt_stats
# ===========================================================================
def run_cutadapt_stats() -> None:
    """
    Parse cutadapt demultiplexing logs and barcodes.csv to build a comprehensive
    JSON & CSV dataset of demultiplexing performance, sample balance, and
    adapter trimming statistics inside METRICS_DIR.
    """
    os.makedirs(METRICS_DIR, exist_ok=True)

    log_candidates = [
        os.path.join(DEMUX_DIR, "log", f"{config.PROJECT_NAME}_step1_cutadapt.log"),
        os.path.join(DEMUX_DIR, "log", f"{config.PROJECT_NAME}_step1_demux.log"),
    ]
    log_file = None
    for cand in log_candidates:
        if os.path.exists(cand):
            log_file = cand
            break

    if not log_file:
        globbed = glob.glob(os.path.join(DEMUX_DIR, "log", "*.log"))
        if globbed:
            log_file = globbed[0]

    if not log_file:
        print(f"  [ERROR] No cutadapt log file found in {os.path.join(DEMUX_DIR, 'log')}. Run step_1 first.")
        return

    print(f"  Parsing cutadapt log file: {log_file}")
    with open(log_file, "r", encoding="utf-8", errors="replace") as fh:
        log_text = fh.read()

    # ---- 1. Overall Summary -------------------------------------------------
    total_processed = 0
    reads_with_adapters = 0
    reads_discarded = 0
    total_bp_processed = 0
    total_bp_written = 0

    m = re.search(r"Total reads processed:\s+([0-9,]+)", log_text)
    if m: total_processed = int(m.group(1).replace(",", ""))

    m = re.search(r"Reads with adapters:\s+([0-9,]+)", log_text)
    if m: reads_with_adapters = int(m.group(1).replace(",", ""))

    m = re.search(r"Reads discarded as untrimmed:\s+([0-9,]+)", log_text)
    if m: reads_discarded = int(m.group(1).replace(",", ""))

    m = re.search(r"Total basepairs processed:\s+([0-9,]+)\s+bp", log_text)
    if m: total_bp_processed = int(m.group(1).replace(",", ""))

    m = re.search(r"Total written \(filtered\):\s+([0-9,]+)\s+bp", log_text)
    if m: total_bp_written = int(m.group(1).replace(",", ""))

    pct_demux = round((reads_with_adapters / total_processed * 100.0), 2) if total_processed > 0 else 0.0
    pct_discarded = round((reads_discarded / total_processed * 100.0), 2) if total_processed > 0 else 0.0

    # ---- 2. Barcode/Sample metadata & per-sample counts ---------------------
    sample_barcodes = {}
    if os.path.exists(config.BARCODES):
        with open(config.BARCODES, newline="", encoding="utf-8-sig") as fh:
            reader = csv.reader(fh)
            if getattr(config, "HAS_HEADER", True):
                next(reader, None)
            for row in reader:
                if not row or len(row) < max(config.SAMPLE_NAME_COL, config.BARCODE_SEQUENCE_COL):
                    continue
                sname = row[config.SAMPLE_NAME_COL - 1].strip()
                bseq  = row[config.BARCODE_SEQUENCE_COL - 1].strip()
                bid   = row[config.BARCODE_ID_COL - 1].strip() if len(row) >= config.BARCODE_ID_COL else sname
                sample_barcodes[sname] = {"sequence": bseq, "id": bid}

    sample_counts = {}
    adapter_sections = re.findall(r"===\s+Adapter\s+([A-Za-z0-9_-]+)\s+===.*?Trimmed:\s+([0-9,]+)\s+times", log_text, re.DOTALL)
    for aname, cnt_str in adapter_sections:
        sample_counts[aname] = int(cnt_str.replace(",", ""))

    per_sample_lines = re.findall(r"\s+([A-Za-z0-9_-]+)\.fastq\.gz\s+([0-9,]+)\s+reads", log_text)
    for sname, cnt_str in per_sample_lines:
        sample_counts[sname] = int(cnt_str.replace(",", ""))

    sample_table = []
    total_trimmed_sum = sum(sample_counts.values())
    all_names = sorted(set(sample_barcodes.keys()) | set(sample_counts.keys()))

    for sname in all_names:
        cnt = sample_counts.get(sname, 0)
        pct_of_trimmed = round((cnt / total_trimmed_sum * 100.0), 2) if total_trimmed_sum > 0 else 0.0
        pct_of_total   = round((cnt / total_processed * 100.0), 3) if total_processed > 0 else 0.0
        binfo = sample_barcodes.get(sname, {"sequence": "N/A", "id": "N/A"})
        sample_table.append({
            "sample_name": sname,
            "barcode_id": binfo["id"],
            "barcode_sequence": binfo["sequence"],
            "read_count": cnt,
            "pct_of_demuxed_reads": pct_of_trimmed,
            "pct_of_total_raw_reads": pct_of_total,
        })

    sample_table.sort(key=lambda x: x["read_count"], reverse=True)

    # ---- 3. Uniformity & Library Balance ------------------------------------
    counts_list = [row["read_count"] for row in sample_table]
    uniformity  = calc_uniformity_metrics(counts_list)

    # ---- 4. Save metrics ----------------------------------------------------
    metrics_data = {
        "project_name": config.PROJECT_NAME,
        "cutadapt_log_path": log_file,
        "overall_summary": {
            "total_reads_processed": total_processed,
            "reads_demultiplexed": reads_with_adapters,
            "pct_demultiplexed": pct_demux,
            "reads_discarded_untrimmed": reads_discarded,
            "pct_discarded": pct_discarded,
            "total_basepairs_processed": total_bp_processed,
            "total_basepairs_written": total_bp_written,
        },
        "sample_uniformity": uniformity,
        "samples": sample_table,
    }

    json_path = os.path.join(METRICS_DIR, f"{config.PROJECT_NAME}_step1_demux_metrics.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(metrics_data, fh, indent=2)

    csv_path = os.path.join(METRICS_DIR, f"{config.PROJECT_NAME}_step1_sample_counts.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_name", "barcode_id", "barcode_sequence", "read_count", "pct_of_demuxed_reads", "pct_of_total_raw_reads"])
        writer.writeheader()
        writer.writerows(sample_table)

    print(f"  ✓ Saved demultiplexing JSON metrics: {json_path}")
    print(f"  ✓ Saved sample counts CSV table  : {csv_path}")
    print(f"    Total Processed: {total_processed:,} | Demultiplexed: {reads_with_adapters:,} ({pct_demux}%) | Discarded: {reads_discarded:,} ({pct_discarded}%)")
    print(f"    Uniformity Gini: {uniformity['gini']} | CV: {uniformity['cv']} | Min/Max fold: {uniformity['min_max_fold']}x")


# ===========================================================================
# CLI & Main
# ===========================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="step_1qc.py — Step 1 Demultiplexing QC job submission and metrics extraction.")
    parser.add_argument("--only", dest="only_step", choices=[s[0] for s in QC_STEPS], help="Run only this step.")
    return parser.parse_args()


def main() -> None:
    args  = parse_args()
    steps = QC_STEPS[:]
    if args.only_step:
        steps = [s for s in steps if s[0] == args.only_step]

    print("=" * 60)
    print(f"  {config.PROJECT_NAME} — step_1qc.py")
    print(f"  Steps: {', '.join(s[0] for s in steps)}")
    print("=" * 60)

    runners = {
        "fastqc":         run_fastqc,
        "cutadapt_stats": run_cutadapt_stats,
    }

    for step_key, desc in steps:
        print(f"\n{'='*60}")
        print(f"[{step_key}] {desc}")
        print("=" * 60)
        t0 = time.time()
        try:
            runners[step_key]()
        except Exception as exc:
            print(f"  [{step_key}] ERROR: {exc}")
        print(f"[{step_key}] Done in {time.time() - t0:.1f}s")

    print("\n" + "=" * 60)
    print("  step_1qc.py complete.")
    print("  Wait for all FastQC bsub jobs to finish, then run:")
    print("    python3 Scripts/step_1qc_report.py")
    print("  Monitor with: bjobs -u $USER")
    print("=" * 60)


if __name__ == "__main__":
    main()
