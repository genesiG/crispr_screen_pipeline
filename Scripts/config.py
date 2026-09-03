#!/usr/bin/env python3

"""
config.py - Configuration module for the genomic analysis pipeline
          - Change file paths for your specific project
"""

import os
import glob

### SETTINGS (Change as needed)
USERNAME = "<your_username>" # Your HPC username
PROJECT_NAME = "CRISPR_SCREEN_PROJECT" # no spaces allowed
PARENT_FOLDER = "CRISPR_SCREEN_PROJECT"
IS_PAIRED_END = False
EXTENSIONS = ('.fastq.gz', '.fastq')
ALLOW_MISMATCHES = True
BARCODE_MISMATCHES = 0
TARGET_MISMATCHES = 1    # 1 mismatch is standard for 20 bp CRISPR guides (2 is too permissive)
# Search for target sequences before and after barcode sequences?
TARGET_BEFORE_BARCODE = False
###

### sgRNA Library Settings (Change as needed)
GENE_SYMBOL_COL = 1 # column index containing gene symbols
ID_COL = 2 # column index containing gRNA IDs (optional, None if not present)
TARGET_SEQUENCE_COL = 3 # column index containing the sgRNA target sequences
###

### Barcode Settings (Change as needed)
HAS_HEADER = True # are there column names at the file?
SAMPLE_NAME_COL = 1
BARCODE_ID_COL = 2
BARCODE_SEQUENCE_COL = 3
###

### ANALYSIS SETTINGS (Change as needed)

# Set the list of all possible TREATMENTS vs CONTROLS comparisons
# Pattern match -- should match the sample names in the count table from step 1
# (example names below -- replace with your own cell-line/condition labels)
CONTROLS = ["CellLineA_baseline", "CellLineA_low_signal_1st", "CellLineA_mid_signal_1st", "CellLineA_high_signal_1st",
            "CellLineB_baseline", "CellLineB_low_signal_1st", "CellLineB_mid_signal_1st", "CellLineB_high_signal_1st",
            "CellLineC_baseline", "CellLineC_low_signal_1st", "CellLineC_mid_signal_1st", "CellLineC_high_signal_1st",
            "CellLineD_baseline", "CellLineD_low_signal_1st", "CellLineD_mid_signal_1st", "CellLineD_high_signal_1st",
            "CellLineA_baseline", "CellLineA_low_signal_2nd", "CellLineA_mid_signal_2nd", "CellLineA_high_signal_2nd",
            "CellLineB_baseline", "CellLineB_low_signal_2nd", "CellLineB_mid_signal_2nd", "CellLineB_high_signal_2nd",
            "CellLineC_baseline", "CellLineC_low_signal_2nd", "CellLineC_mid_signal_2nd", "CellLineC_high_signal_2nd",
            "CellLineD_baseline", "CellLineD_low_signal_2nd", "CellLineD_mid_signal_2nd", "CellLineD_high_signal_2nd"]
TREATMENTS = ["CellLineA_low_signal_1st", "CellLineA_mid_signal_1st", "CellLineA_high_signal_1st",
              "CellLineB_low_signal_1st", "CellLineB_mid_signal_1st", "CellLineB_high_signal_1st",
              "CellLineC_low_signal_1st", "CellLineC_mid_signal_1st", "CellLineC_high_signal_1st",
              "CellLineD_low_signal_1st", "CellLineD_mid_signal_1st", "CellLineD_high_signal_1st",
              "CellLineA_low_signal_2nd", "CellLineA_mid_signal_2nd", "CellLineA_high_signal_2nd",
              "CellLineB_low_signal_2nd", "CellLineB_mid_signal_2nd", "CellLineB_high_signal_2nd",
              "CellLineC_low_signal_2nd", "CellLineC_mid_signal_2nd", "CellLineC_high_signal_2nd",
              "CellLineD_low_signal_2nd", "CellLineD_mid_signal_2nd", "CellLineD_high_signal_2nd"]
EXCLUDE_SAMPLES = ["xyz"]  # String patterns to search for samples to exclude (optional)



### Work directory set up
HOMEDIR = os.path.join("/home", USERNAME)
WORKDIR = os.path.join(HOMEDIR, PARENT_FOLDER)
CODEDIR = os.path.join(WORKDIR, "Scripts")
ORIGINAL_DATA = os.path.join(WORKDIR, "Original_Data")
METADATA = os.path.join(WORKDIR, "Metadata")
ANALYSIS_DATA = os.path.join(WORKDIR, "Analysis_Data")
IMPORTABLE_DATA = os.path.join(WORKDIR, "Importable_Data")
QC_DIR = os.path.join(ANALYSIS_DATA, "qc")
MULTIQCDIR = os.path.join(QC_DIR, "multiqc")
###

### REFERENCE FILES (Change as needed)
# Path to folder with raw sequencing files
DATADIR = os.path.join(ORIGINAL_DATA, "fastq_files")
# Set path for fastq file
FASTQ_PATH = os.path.join(DATADIR, "your_sequencing_run.fastq.gz")
# Path to barcodes file
BARCODES = os.path.join(METADATA, "barcodes.csv")
# Path to sgRNA library file
LIBRARY = glob.glob(os.path.join(METADATA, "grna_library_pool*.csv"))
# Path to NTC gRNAs
NTC_FILE = os.path.join(METADATA, "ntc_sgRNA_library.csv")
# Output FASTQC reports
QCDIR1 = os.path.join(ANALYSIS_DATA, "fastqc1")
# Output read counts data
COUNTDIR = os.path.join(ANALYSIS_DATA, "counts")
###
