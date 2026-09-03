#!/usr/bin/env python3
"""
step_2qc.py
===========
Guide representation and QC plotting — submits an LSF batch job that runs the
`guide_representation.R` helper script over the MAGeCK count tables produced
by step_2_count.py.

Tool used : Rscript Scripts/utils/guide_representation.R
LSF       : 1 job, 4 cores, 2 h wall time

Prerequisites:
    Run step_2_count.py and wait for the count tables to be generated before
    running this step.

Skip-if-exists:
    If plots already exist in Analysis_Data/counts/plots/, prints status.

Inputs (from Scripts/config.py):
    config.COUNTDIR/*.count.txt  — MAGeCK count tables from step 2

Outputs:
    Analysis_Data/counts/plots/*_guide_representation.png
    Analysis_Data/counts/plots/*_ecdf.png
    Analysis_Data/counts/plots/*_lorenz.png
    Analysis_Data/counts/plots/*_dropout_summary.csv
    Analysis_Data/counts/plots/*_density.png
    Analysis_Data/counts/plots/*_controls_ecdf.png (if control guides exist)
    Analysis_Data/counts/log/<PROJECT>_step2qc.log
    Analysis_Data/counts/log/<PROJECT>_step2qc.error
"""

import glob
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


def main() -> None:
    print("=" * 60)
    print(f"Project : {config.PROJECT_NAME}")
    print(f"Step    : 2 QC (Guide Representation & Lorenz/Gini/ECDF Plots)")
    print("=" * 60)

    # ---- Check inputs & setup paths -----------------------------------------
    count_files = glob.glob(os.path.join(config.COUNTDIR, "*.count.txt"))
    if not count_files:
        print(
            f"\n[WARNING] No *.count.txt files found in {config.COUNTDIR}.\n"
            "If step_2_count.py is currently running on LSF, you can still submit this job\n"
            "or wait for step_2_count.py to finish before running step_2qc.py.\n"
        )
    else:
        print(f"Found {len(count_files)} count table(s) in {config.COUNTDIR}:")
        for cf in sorted(count_files):
            print(f"  ✓ {os.path.basename(cf)}")

    helper_script = os.path.join(os.path.dirname(__file__), "utils", "guide_representation.R")
    if not os.path.exists(helper_script):
        print(f"[ERROR] Helper script not found: {helper_script}", file=sys.stderr)
        sys.exit(1)

    log_dir = os.path.join(config.COUNTDIR, "log")
    os.makedirs(log_dir, exist_ok=True)
    plot_dir = os.path.join(config.COUNTDIR, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    log_out = os.path.join(log_dir, f"{config.PROJECT_NAME}_step2qc.log")
    log_err = os.path.join(log_dir, f"{config.PROJECT_NAME}_step2qc.error")

    job_script = f"""#!/bin/bash
#BSUB -J "step2_qc_{config.PROJECT_NAME}"
#BSUB -q normal
#BSUB -o "{log_out}"
#BSUB -e "{log_err}"
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -W 2:00

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate crisprscreen

set -euo pipefail

echo "=== Step 2 QC: Guide Representation Analysis ==="
echo "Input Dir : {config.COUNTDIR}"
echo "Output Dir: {plot_dir}"
echo "Started   : $(date)"

# Pass all *.count.txt files in COUNTDIR directly to guide_representation.R
Rscript "{helper_script}" {config.COUNTDIR}/*.count.txt 2>&1 | tee -a "{log_out}"

echo "Finished  : $(date)"
echo "=== Generated QC Plots & Tables ==="
ls -lh "{plot_dir}"
"""

    print("\nSubmitting guide representation QC job (Rscript)...")
    _bsub(job_script, label="step2_qc")
    print(
        "\nMonitor with:  bjobs -u $USER\n"
        f"Plots will be written to:  {plot_dir}"
    )


if __name__ == "__main__":
    main()
