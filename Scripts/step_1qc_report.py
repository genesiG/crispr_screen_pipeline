#!/usr/bin/env python3
"""
step_1qc_report.py
==================
Step 1 Demultiplexing QC Report Generator — runs after all FastQC jobs from
step_1qc.py have finished.

Steps (sequential — each step has its own subdirectory inside config.QC_DIR):
  multiqc     → MultiQC aggregation of FastQC reports (bsub+wait) -> qc/multiqc
  demux_plots → Generate demultiplexing read count, fate & uniformity plots -> qc/demux_plots
  report      → Generate combined self-contained HTML QC report -> qc/report

Usage:
    python3 Scripts/step_1qc_report.py                 # run all steps
    python3 Scripts/step_1qc_report.py --only report   # regenerate HTML report only
"""

import argparse
import datetime
import json
import math
import os
import re
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# Import pipeline config (Scripts/ directory)
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
import config

REPORT_STEPS = [
    ("multiqc",     "MultiQC aggregation of demultiplexed FastQC reports (bsub+wait) -> qc/multiqc"),
    ("demux_plots", "Generate demultiplexing read count, fate & uniformity plots -> qc/demux_plots"),
    ("report",      "Generate combined self-contained HTML QC report -> qc/report"),
]

FASTQC_DIR      = os.path.join(config.QC_DIR, "fastqc")
METRICS_DIR     = os.path.join(config.QC_DIR, "cutadapt_stats")
PLOTS_DIR       = os.path.join(config.QC_DIR, "demux_plots")
MULTIQC_DIR     = config.MULTIQCDIR
REPORT_DIR      = os.path.join(config.QC_DIR, "report")
REPORT_PATH     = os.path.join(REPORT_DIR, f"{config.PROJECT_NAME}_step1_qc_report.html")
METRICS_JSON    = os.path.join(METRICS_DIR, f"{config.PROJECT_NAME}_step1_demux_metrics.json")


# ===========================================================================
# Shared bsub polling helper
# ===========================================================================
def _submit_bsub_and_wait(batch_cmd: str, log_dir: str, step_key: str,
                          poll_interval: int = 15, timeout_minutes: int = 60) -> bool:
    os.makedirs(log_dir, exist_ok=True)
    batch_file = os.path.join(log_dir, f"{step_key}.batch")
    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    result = subprocess.run(["bsub"], input=batch_cmd, text=True, capture_output=True, check=False)
    stdout = result.stdout.strip()
    print(f"  bsub output: {stdout}")

    m = re.search(r"Job <(\d+)>", stdout)
    if not m:
        print("  [WARNING] Could not parse job ID from bsub output. Proceeding without polling.")
        return True

    job_id = m.group(1)
    print(f"  Job submitted: {job_id}. Polling every {poll_interval}s...")

    deadline = time.time() + timeout_minutes * 60
    while time.time() < deadline:
        time.sleep(poll_interval)
        bjobs = subprocess.run(["bjobs", "-noheader", job_id], text=True, capture_output=True, check=False)
        status_line = bjobs.stdout.strip()
        if not status_line:
            bhist = subprocess.run(["bhist", "-noheader", "-d", job_id], text=True, capture_output=True, check=False)
            if "DONE" in bhist.stdout or bhist.stdout.strip() == "":
                print(f"  Job {job_id} DONE.")
                return True
            if "EXIT" in bhist.stdout:
                raise RuntimeError(f"Job {job_id} ({step_key}) exited with failure.")
            print(f"  Job {job_id} no longer in queue (assumed DONE).")
            return True

        fields = status_line.split()
        if len(fields) >= 3:
            state = fields[2]
            if state == "DONE":
                print(f"  Job {job_id} DONE.")
                return True
            if state == "EXIT":
                raise RuntimeError(f"Job {job_id} ({step_key}) exited with failure.")
            print(f"  Job {job_id} status: {state} ...")

    print(f"  [WARNING] Timed out waiting for job {job_id} after {timeout_minutes} min.")
    return False


# ===========================================================================
# multiqc — aggregate demultiplexed FastQC reports inside qc/multiqc
# ===========================================================================
def run_multiqc() -> None:
    """Run MultiQC over FASTQC_DIR to aggregate demultiplexed sample FastQC reports into MULTIQC_DIR."""
    log_dir = os.path.join(MULTIQC_DIR, "log")
    os.makedirs(MULTIQC_DIR, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    job_name     = f"{config.PROJECT_NAME}_step1_multiqc"
    log_out      = os.path.join(log_dir, "multiqc.log")
    log_err      = os.path.join(log_dir, "multiqc.error")
    multiqc_html = f"{config.PROJECT_NAME}_step1_multiqc_report.html"

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

echo "=== Running MultiQC on {FASTQC_DIR} ==="
multiqc "{FASTQC_DIR}" --outdir "{MULTIQC_DIR}" --filename "{multiqc_html}" --force --verbose
echo "=== MultiQC finished successfully ==="
"""
    _submit_bsub_and_wait(batch_cmd, log_dir, "multiqc")


# ===========================================================================
# demux_plots — generate HTML/SVG plots of demultiplexing metrics inside qc/demux_plots
# ===========================================================================
def run_demux_plots() -> None:
    """Generate high-quality HTML/SVG plots for read counts, fate, and uniformity inside PLOTS_DIR."""
    os.makedirs(PLOTS_DIR, exist_ok=True)

    if not os.path.exists(METRICS_JSON):
        print(f"  [INFO] {METRICS_JSON} not found. Running cutadapt stats extraction...")
        import step_1qc
        step_1qc.run_cutadapt_stats()

    if not os.path.exists(METRICS_JSON):
        print(f"  [WARNING] Could not find or create {METRICS_JSON}. Skipping demux_plots.")
        return

    with open(METRICS_JSON, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    samples = data.get("samples", [])
    summary = data.get("overall_summary", {})

    if not samples:
        print("  [WARNING] No samples found in metrics JSON.")
        return

    # --- 1. Barplot of read counts per sample (SVG) ---
    max_count = max([s["read_count"] for s in samples] or [1])
    svg_bars = []
    bar_height = 36
    top_margin = 60
    chart_height = top_margin + len(samples) * bar_height + 40
    chart_width = 800
    bar_max_w = 480

    for i, s in enumerate(samples):
        y = top_margin + i * bar_height
        w = int((s["read_count"] / max_count) * bar_max_w) if max_count > 0 else 0
        color = "#88419d" if i % 2 == 0 else "#6e016b"
        svg_bars.append(f"""
        <g transform="translate(0, {y})">
          <text x="200" y="18" fill="#e2e8f0" font-family="sans-serif" font-size="13" text-anchor="end" font-weight="600">{s['sample_name']}</text>
          <rect x="210" y="3" width="{w}" height="22" fill="{color}" rx="4" />
          <text x="{218 + w}" y="18" fill="#a0aec0" font-family="sans-serif" font-size="12">{s['read_count']:,} ({s['pct_of_demuxed_reads']}%)</text>
        </g>""")

    svg_read_counts = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {chart_width} {chart_height}" width="100%" height="{chart_height}px" style="background:#13161f; border-radius:8px;">
      <text x="20" y="35" fill="#e2e8f0" font-family="sans-serif" font-size="16" font-weight="700">Per-Sample Demultiplexed Read Counts</text>
      {"".join(svg_bars)}
    </svg>"""

    out_svg_bars = os.path.join(PLOTS_DIR, f"{config.PROJECT_NAME}_demux_read_counts.svg")
    with open(out_svg_bars, "w", encoding="utf-8") as fh:
        fh.write(svg_read_counts)

    # --- 2. Read Fate Breakdown Pie Chart (SVG) ---
    demuxed = summary.get("reads_demultiplexed", 0)
    discarded = summary.get("reads_discarded_untrimmed", 0)
    total = summary.get("total_reads_processed", demuxed + discarded)
    pct_demux = summary.get("pct_demultiplexed", 0.0)
    pct_disc = summary.get("pct_discarded", 0.0)

    angle = (demuxed / total * 360.0) if total > 0 else 0
    rad = math.radians(angle)
    x_end = 200 + 100 * math.sin(rad)
    y_end = 200 - 100 * math.cos(rad)
    large_arc = 1 if angle > 180 else 0

    donut_svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 700 400" width="100%" height="400px" style="background:#13161f; border-radius:8px;">
      <text x="20" y="35" fill="#e2e8f0" font-family="sans-serif" font-size="16" font-weight="700">Read Fate Breakdown (cutadapt demultiplexing)</text>
      <!-- Background circle (Discarded) -->
      <circle cx="200" cy="200" r="100" fill="#374151" />
      <!-- Demultiplexed slice -->
      <path d="M 200 200 L 200 100 A 100 100 0 {large_arc} 1 {x_end:.1f} {y_end:.1f} Z" fill="#88419d" />
      <!-- Center hole for donut -->
      <circle cx="200" cy="200" r="55" fill="#13161f" />
      <!-- Legend & Stats -->
      <g transform="translate(360, 150)">
        <rect x="0" y="0" width="16" height="16" fill="#88419d" rx="3" />
        <text x="26" y="13" fill="#e2e8f0" font-family="sans-serif" font-size="14" font-weight="600">Demultiplexed Reads: {demuxed:,} ({pct_demux}%)</text>
        <rect x="0" y="36" width="16" height="16" fill="#374151" rx="3" />
        <text x="26" y="49" fill="#e2e8f0" font-family="sans-serif" font-size="14" font-weight="600">Discarded / Untrimmed: {discarded:,} ({pct_disc}%)</text>
        <text x="0" y="90" fill="#a0aec0" font-family="sans-serif" font-size="13">Total Processed: {total:,} reads</text>
      </g>
    </svg>"""

    out_svg_donut = os.path.join(PLOTS_DIR, f"{config.PROJECT_NAME}_demux_read_fate.svg")
    with open(out_svg_donut, "w", encoding="utf-8") as fh:
        fh.write(donut_svg)

    print(f"  ✓ Generated SVG plots in {PLOTS_DIR}")


# ===========================================================================
# report — Build combined HTML QC report inside qc/report
# ===========================================================================
_HTML_STYLE = """
<style>
  :root {
    --bg:      #0f1117;
    --surface: #1a1d27;
    --border:  #2d3148;
    --accent:  #88419d;
    --text:    #e2e8f0;
    --muted:   #6b7280;
    --warn:    #f59e0b;
    --font:    'Inter', 'Segoe UI', system-ui, sans-serif;
  }
  @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--font);
    font-size: 15px;
    line-height: 1.6;
    padding: 2rem;
    max-width: 1400px;
    margin: 0 auto;
  }
  header { border-bottom: 2px solid var(--accent); padding-bottom: 1.5rem; margin-bottom: 2.5rem; }
  header h1 { font-size: 2rem; font-weight: 700; color: var(--text); letter-spacing: -0.02em; }
  header p { color: var(--muted); margin-top: 0.4rem; }
  nav { background: var(--surface); border: 1px solid var(--border); border-radius: 8px; padding: 1rem 1.5rem; margin-bottom: 2.5rem; }
  nav h3 { font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.1em; color: var(--muted); margin-bottom: 0.6rem; }
  nav ol { padding-left: 1.4rem; }
  nav ol li { margin-bottom: 0.25rem; }
  nav ol li a { color: var(--accent); text-decoration: none; }
  nav ol li a:hover { text-decoration: underline; }
  .qc-section { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 1.8rem 2rem; margin-bottom: 2.5rem; }
  .qc-section h2 { font-size: 1.25rem; font-weight: 700; color: var(--accent); padding-bottom: 0.6rem; border-bottom: 1px solid var(--border); margin-bottom: 1.2rem; }
  .kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1.2rem; margin-bottom: 2rem; }
  .kpi-card { background: #13161f; border: 1px solid var(--border); border-radius: 8px; padding: 1.2rem; text-align: center; }
  .kpi-value { font-size: 1.6rem; font-weight: 700; color: #fff; margin: 0.4rem 0; }
  .kpi-label { font-size: 0.85rem; color: var(--muted); text-transform: uppercase; letter-spacing: 0.05em; }
  table { width: 100%; border-collapse: collapse; margin: 1.5rem 0; font-size: 0.92rem; }
  th, td { padding: 0.75rem 1rem; border-bottom: 1px solid var(--border); text-align: left; }
  th { background: #13161f; color: var(--text); font-weight: 600; }
  tr:hover { background: rgba(136, 65, 157, 0.1); }
  figure { margin: 1.2rem 0; padding: 0.5rem; background: #13161f; border-radius: 6px; border: 1px solid var(--border); }
  footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid var(--border); color: var(--muted); font-size: 0.8rem; text-align: center; }
</style>
"""


def _svg_inline(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def run_report() -> None:
    """Compile all demultiplexing stats, tables, SVG plots, and MultiQC into one HTML document inside REPORT_DIR."""
    os.makedirs(REPORT_DIR, exist_ok=True)
    print(f"  Generating combined HTML QC report → {REPORT_PATH}")

    if not os.path.exists(METRICS_JSON):
        print(f"  [WARNING] {METRICS_JSON} not found. Running cutadapt_stats first...")
        import step_1qc
        step_1qc.run_cutadapt_stats()

    data = {}
    if os.path.exists(METRICS_JSON):
        with open(METRICS_JSON, "r", encoding="utf-8") as fh:
            data = json.load(fh)

    samples = data.get("samples", [])
    summary = data.get("overall_summary", {})
    uniformity = data.get("sample_uniformity", {})

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # --- Section 1: Executive KPI Summary ---
    kpi_html = f"""<div class="kpi-grid">
      <div class="kpi-card">
        <div class="kpi-label">Total Processed Reads</div>
        <div class="kpi-value">{summary.get('total_reads_processed', 0):,}</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Reads Demultiplexed</div>
        <div class="kpi-value" style="color:#10b981;">{summary.get('reads_demultiplexed', 0):,}</div>
        <div style="font-size:0.8rem; color:var(--muted);">{summary.get('pct_demultiplexed', 0)}% of total</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Sample Uniformity (Gini)</div>
        <div class="kpi-value" style="color:{'#10b981' if uniformity.get('gini', 1) < 0.2 else '#f59e0b'};">{uniformity.get('gini', 'N/A')}</div>
        <div style="font-size:0.8rem; color:var(--muted);">CV: {uniformity.get('cv', 'N/A')}</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Min/Max Pool Ratio</div>
        <div class="kpi-value">{uniformity.get('min_max_fold', 'N/A')}x</div>
        <div style="font-size:0.8rem; color:var(--muted);">Median: {uniformity.get('median', 0):,} reads</div>
      </div>
    </div>"""

    # --- Section 2: Cutadapt Diagnostic Summary Table ---
    diag_table = f"""<table>
      <thead>
        <tr><th>Metric</th><th>Count / Value</th><th>Percentage</th><th>Description</th></tr>
      </thead>
      <tbody>
        <tr><td>Total Reads Processed</td><td>{summary.get('total_reads_processed', 0):,}</td><td>100.0%</td><td>Total raw reads in input FASTQ ({os.path.basename(config.FASTQ_PATH)})</td></tr>
        <tr><td>Reads Demultiplexed (Passing Filters)</td><td style="color:#10b981; font-weight:600;">{summary.get('reads_demultiplexed', 0):,}</td><td>{summary.get('pct_demultiplexed', 0)}%</td><td>Reads with valid exact barcode matches written to per-sample FASTQs</td></tr>
        <tr><td>Discarded as Untrimmed</td><td style="color:#f59e0b;">{summary.get('reads_discarded_untrimmed', 0):,}</td><td>{summary.get('pct_discarded', 0)}%</td><td>Reads lacking barcode match (PhiX spike-in, sequencing errors, or off-target reads)</td></tr>
        <tr><td>Total Basepairs Processed</td><td>{summary.get('total_basepairs_processed', 0):,} bp</td><td>-</td><td>Total nucleotide volume processed</td></tr>
        <tr><td>Total Basepairs Written</td><td>{summary.get('total_basepairs_written', 0):,} bp</td><td>-</td><td>Basepair volume retained across demultiplexed samples</td></tr>
      </tbody>
    </table>"""

    # --- Section 3: Per-Sample Read Counts Table ---
    rows_html = []
    for s in samples:
        rows_html.append(f"""<tr>
          <td style="font-weight:600; color:#fff;">{s['sample_name']}</td>
          <td><code>{s['barcode_id']}</code></td>
          <td><code style="color:#a855f7;">{s['barcode_sequence']}</code></td>
          <td style="text-align:right; font-weight:600;">{s['read_count']:,}</td>
          <td style="text-align:right;">{s['pct_of_demuxed_reads']}%</td>
          <td style="text-align:right;">{s['pct_of_total_raw_reads']}%</td>
        </tr>""")

    sample_table = f"""<table>
      <thead>
        <tr><th>Sample / Library Pool</th><th>Barcode ID</th><th>Barcode Sequence</th><th style="text-align:right;">Read Count</th><th style="text-align:right;">% of Demultiplexed</th><th style="text-align:right;">% of Total Raw</th></tr>
      </thead>
      <tbody>
        {"".join(rows_html)}
      </tbody>
    </table>"""

    # --- Section 4: Plots ---
    svg_bars_path  = os.path.join(PLOTS_DIR, f"{config.PROJECT_NAME}_demux_read_counts.svg")
    svg_donut_path = os.path.join(PLOTS_DIR, f"{config.PROJECT_NAME}_demux_read_fate.svg")

    plots_html = ""
    if os.path.exists(svg_bars_path):
        plots_html += f'<figure class="svg-figure">{_svg_inline(svg_bars_path)}</figure>'
    if os.path.exists(svg_donut_path):
        plots_html += f'<figure class="svg-figure">{_svg_inline(svg_donut_path)}</figure>'

    # --- Section 5: MultiQC Embed ---
    multiqc_file = os.path.join(MULTIQC_DIR, f"{config.PROJECT_NAME}_step1_multiqc_report.html")
    multiqc_html = ""
    if os.path.exists(multiqc_file):
        try:
            with open(multiqc_file, encoding="utf-8") as fh:
                raw_mqc = fh.read()
            srcdoc = raw_mqc.replace("&", "&amp;").replace('"', "&quot;")
            multiqc_html = f"""<figure>
              <iframe srcdoc="{srcdoc}" width="100%" height="700px" frameborder="0" sandbox="allow-scripts allow-same-origin allow-popups">
                <p><a href="{multiqc_file}">Open MultiQC Report</a></p>
              </iframe>
              <p style="font-size:0.85rem; color:var(--muted); margin-top:0.5rem;"><a href="{multiqc_file}" target="_blank">↗ Open Full Interactive MultiQC Report in New Tab</a></p>
            </figure>"""
        except Exception as e: # pylint: disable=broad-except
            multiqc_html = f'<p class="missing">⚠ Could not embed MultiQC report: {e}</p>'
    else:
        multiqc_html = f'<p class="missing">⚠ MultiQC report not found at <code>{multiqc_file}</code>. Run multiqc step first.</p>'

    # --- Assemble HTML ---
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{config.PROJECT_NAME} — Step 1 Demultiplexing QC Report</title>
  {_HTML_STYLE}
</head>
<body>

<header>
  <h1>🧬 {config.PROJECT_NAME} — Step 1 Demultiplexing QC Report</h1>
  <p>Generated: {now} &nbsp;|&nbsp; Input FASTQ: <code>{os.path.basename(config.FASTQ_PATH)}</code> &nbsp;|&nbsp; Barcodes: <code>{os.path.basename(config.BARCODES)}</code></p>
</header>

<nav>
  <h3>Contents</h3>
  <ol>
    <li><a href="#kpi">Executive Summary & Uniformity KPIs</a></li>
    <li><a href="#fate">Cutadapt Read Fate Diagnostics</a></li>
    <li><a href="#table">Per-Sample Read Counts & Barcode Summary</a></li>
    <li><a href="#plots">Demultiplexing Performance Plots</a></li>
    <li><a href="#multiqc">MultiQC Per-Sample Quality Metrics</a></li>
  </ol>
</nav>

<section class="qc-section" id="kpi">
  <h2>1. Executive Summary & Uniformity KPIs</h2>
  {kpi_html}
</section>

<section class="qc-section" id="fate">
  <h2>2. Cutadapt Read Fate Diagnostics</h2>
  {diag_table}
</section>

<section class="qc-section" id="table">
  <h2>3. Per-Sample Read Counts & Barcode Summary</h2>
  {sample_table}
</section>

<section class="qc-section" id="plots">
  <h2>4. Demultiplexing Performance Plots</h2>
  {plots_html}
</section>

<section class="qc-section" id="multiqc">
  <h2>5. MultiQC Per-Sample Quality Metrics</h2>
  {multiqc_html}
</section>

<footer>
  <p>Generated by <code>step_1qc_report.py</code> — {config.PROJECT_NAME}</p>
</footer>

</body>
</html>
"""

    with open(REPORT_PATH, "w", encoding="utf-8") as fh:
        fh.write(html)

    print(f"  ✓ Combined HTML report written to: {REPORT_PATH}")


# ===========================================================================
# CLI & Main
# ===========================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="step_1qc_report.py — Compile Step 1 Demultiplexing QC Report.")
    parser.add_argument("--only", dest="only_step", choices=[s[0] for s in REPORT_STEPS], help="Run only this step.")
    return parser.parse_args()


def main() -> None:
    args  = parse_args()
    steps = REPORT_STEPS[:]
    if args.only_step:
        steps = [s for s in steps if s[0] == args.only_step]

    print("=" * 60)
    print(f"  {config.PROJECT_NAME} — step_1qc_report.py")
    print(f"  Steps: {', '.join(s[0] for s in steps)}")
    print("=" * 60)

    runners = {
        "multiqc":     run_multiqc,
        "demux_plots": run_demux_plots,
        "report":      run_report,
    }

    for step_key, desc in steps:
        print(f"\n{'='*60}")
        print(f"[{step_key}] {desc}")
        print("=" * 60)
        t0 = time.time()
        try:
            runners[step_key]()
        except Exception as exc: # pylint: disable=broad-except
            print(f"  [{step_key}] ERROR: {exc}")
        print(f"[{step_key}] Done in {time.time() - t0:.1f}s")

    print("\n" + "=" * 60)
    print("  step_1qc_report.py complete.")
    print(f"  HTML report: {REPORT_PATH}")
    print("=" * 60)


if __name__ == "__main__":
    main()
