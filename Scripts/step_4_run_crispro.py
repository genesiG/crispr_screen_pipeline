#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_4_run_crispro.py - Submit CRISPRO analysis as batch jobs
Steps:
1. Process sgrna_summary files to create crispro input files
2. For each processed file, submit CRISPRO job with appropriate parameters
3. Submit batch jobs to cluster using bsub
"""

import os
import time
import glob
import pandas as pd
import config

# Parameters (no defaults)
fastq_files = []                   # List of FASTQ filenames

# Optional parameters with defaults
conditions = []              # Experimental conditions
comparisons = []             # Comparisons (treatment,control)
scoring_method = "ALFC"      # ALFC or DESEQ2
negative_controls = config.NTC_FILE     # File with negative control guides
positive_controls = None     # File with positive control guides
guide_len = 20               # Length of guides
adapter = None               # Adapter 5' of guide
num_processes = 4            # Number of processes
loess_formula = "position"   # Loess regression modeling
color_structure = True       # Download and recolor structures
keep_temp = True             # Keep intermediate files
principal_iso = False        # Use principal isoform. If false, will map to all isoforms
isoform_list = []            # Prefered transcript IDs
def_hits = []                # Score definitions and hits
offtarget_filter = 5         # Off-target score filter

### Helper Functions ###

def process_sgrna_summary_files(target_gene):
    """
    Process sgrna_summary files from mageck to create crispro input files:
    - Load library file
    - For each sgrna_summary file:
        * Merge with library data
        * Filter for EED|NTC targets
        * Save as crispro input file
    Returns list of generated crispro files
    """
    # Set up directories
    workdir = config.CRISPRODIR
    importable_dir = config.IMPORTABLE_DATA
    lib_path = config.LIBRARY

    # Read library file
    lib = pd.read_csv(lib_path)
    lib = lib[["sgRNA sequence", "Target Gene Symbol", "Target Gene ID"]]
    lib = lib.rename(columns={
        "Target Gene Symbol": "sgrna",
        "sgRNA sequence": "sequence",
        "Target Gene ID": "Gene"
    })

    # Find all sgrna_summary files
    summary_files = glob.glob(os.path.join(importable_dir, config.SGRNA_SUMMARY_PATTERN))
    crispro_files = []

    print(f"Found {len(summary_files)} sgrna_summary files to process")

    for file_path in summary_files:
        # Read and merge data
        results = pd.read_csv(file_path, sep='\t')
        merged = pd.merge(lib, results, on="sgrna")

        # Using case-insensitive matching with target gene and 'NTC' patterns
        gene_mask = merged['sgrna'].str.contains(target_gene, case=False, na=False)
        ntc_mask = merged['sgrna'].str.contains(config.NTC_PATTERN, case=False, na=False)
        merged = merged[gene_mask | ntc_mask][["sequence", "LFC"]]
        # merged = merged[["sequence", "LFC"]]

        # Create output filename
        base_name = os.path.basename(file_path)
        out_name = base_name.replace("sgrna_summary.txt", f"crispro_{target_gene}.csv")
        targetdir = os.path.join(workdir, target_gene)
        scoredir = os.path.join(targetdir, "scores")
        out_path = os.path.join(scoredir, out_name)

        os.makedirs(workdir, exist_ok=True)
        os.makedirs(scoredir, exist_ok=True)
        os.makedirs(targetdir, exist_ok=True)

        # Save processed file
        merged.to_csv(out_path, sep=',', index=False)
        crispro_files.append(out_path)
        print(f"Created crispro input file: {out_path}")

    return crispro_files

def create_crispro_command(score_file, target, job_id):
    """
    Create CRISPRO command components for a specific score file
    Returns command parts list and job name
    """
    # Define output directory
    targetdir = os.path.join(config.CRISPRODIR, f"{target['gene']}_{target['pdb']}")
    os.makedirs(targetdir, exist_ok=True)
    results_dir = os.path.join(targetdir, "results")
    os.makedirs(results_dir, exist_ok=True)
    outdir = os.path.join(results_dir, job_id)

    # Parameters
    identifier = target["identifier"]
    pdb = target["pdb"]

    # Build command components
    cmd_parts = ["crispro"]

    # Required parameters
    if not identifier:
        raise ValueError("Identifier must be provided for CRISPRO analysis")
    cmd_parts.extend(["-a", config.ANNOFILE])
    cmd_parts.extend(["-i", identifier])

    # FASTQ files
    if fastq_files:
        cmd_parts.extend(["-f"] + fastq_files)

    # Optional parameters
    if conditions:
        cmd_parts.extend(["--conditions"] + conditions)

    if comparisons:
        cmd_parts.extend(["--comparisons"] + [f"{t},{c}" for t,c in comparisons])

    cmd_parts.extend(["--scoring-method", scoring_method])

    if negative_controls:
        cmd_parts.extend(["--negative-controls", negative_controls])

    if positive_controls:
        cmd_parts.extend(["--positive-controls", positive_controls])

    cmd_parts.extend(["--guide-len", str(guide_len)])

    if adapter:
        cmd_parts.extend(["--adapter", adapter])

    # Add score file parameter
    cmd_parts.extend(["--scores", score_file])

    cmd_parts.extend(["-o", outdir])
    cmd_parts.extend(["--num-processes", str(num_processes)])
    cmd_parts.extend(["--loess-formula", loess_formula])

    if color_structure:
        cmd_parts.append("--color-structure")

    if keep_temp:
        cmd_parts.append("--keep-temp")

    if principal_iso:
        cmd_parts.append("--principal-iso")

    if isoform_list:
        cmd_parts.extend(["--isoform-list"] + isoform_list)

    if pdb and identifier:
        cmd_parts.extend(["--extra-pdb"] + [f"{pdb},{identifier}"])

    if def_hits:
        cmd_parts.extend(["--def-hits"] + def_hits)

    cmd_parts.extend(["--offtarget-filter", str(offtarget_filter)])

    # Create job name from score file
    base_name = os.path.basename(score_file)
    # job_name = f"crispro_{target['gene']}_{base_name.replace('.crispro.txt', '')}"
    job_name = f"crispro_{target['gene']}_{base_name.split('.')[0]}"

    return cmd_parts, job_name, outdir

def submit_crispro_job(cmd_parts, job_name, outdir):
    """
    Submit CRISPRO job as a batch job
    """
    log_dir = os.path.join(config.CRISPRODIR, "logs")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    # Build batch script content
    batch_content = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -o {log_dir}/{job_name}.log
#BSUB -e {log_dir}/{job_name}.error
#BSUB -q normal
#BSUB -n 8,32
#BSUB -M 64400
#BSUB -R "rusage [mem=64400] span[hosts=1]"

# Load required modules
conda activate crispro

# Execute CRISPRO
{' '.join(cmd_parts)}
"""

    # Write batch script
    with open(batch_file, 'w', encoding='utf-8') as f:
        f.write(batch_content)

    # Submit job
    print(f"Submitting CRISPRO job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

### Main Workflow ###
def main():
    """Main function to submit job"""

    # Process each target
    for target in config.TARGETS:
        print(f"\nProcessing target: {target['gene']}")

        # Step 1: Create input files
        crispro_files = process_sgrna_summary_files(target["gene"])

        if not crispro_files:
            print(f"No files created for {target['gene']}. Skipping.")
            continue

        # Step 2: Submit jobs
        print("\nSubmitting CRISPRO jobs for processed files:")
        for idx, score_file in enumerate(crispro_files):
            print(f"\nProcessing file {idx+1}/{len(crispro_files)}: {score_file}")
            basename = os.path.basename(score_file)
            # job_id = basename.replace(".crispro.txt", "")
            job_id = basename.replace(f".crispro_{target['gene']}.csv", "")
            cmd_parts, job_name, outdir = create_crispro_command(
                score_file, target, job_id
            )
            submit_crispro_job(cmd_parts, job_name, outdir)

    # # Step 1: Process sgrna_summary files
    # print("Processing sgrna_summary files...")
    # crispro_files = process_sgrna_summary_files()

    # if not crispro_files:
    #     print("No crispro files created. Exiting.")
    #     return

    # print("\nSubmitting CRISPRO jobs for processed files:")

    # # Step 2: Submit CRISPRO job for each processed file
    # for idx, score_file in enumerate(crispro_files):
    #     print(f"\nProcessing file {idx+1}/{len(crispro_files)}: {score_file}")
    #     basename = os.path.basename(score_file)
    #     job_id = basename.replace(".crispro.txt", "")

    #     # Create command for this score file
    #     cmd_parts, job_name, outdir = create_crispro_command(score_file, job_id)

    #     # Submit job
    #     submit_crispro_job(cmd_parts, job_name, outdir)

    print("\nAll CRISPRO jobs submitted successfully")

if __name__ == "__main__":
    main()
