# SKILL: skill_plot.md
# PURPOSE: Generate, save, and interpret QC and summary figures for the pipeline.

## When to load this skill
Load when:
- A step produces quantitative output (alignment rates, read counts, QC metrics)
- You need to visualize a summary table or distribution
- The user asks you to "plot", "visualize", or "check QC"
- Alignment rates from step_2 need to be evaluated

---

## General Rules
- Always modify one of the optional scripts (e.g. `step_2b_*`, `step_4b_*`) or create a new dedicated script that fits seamlesly in the pipeline
    - e.g. `Scripts/step_3b_bam_qc.py`    
- All figures go in `*/plots/` INSIDE their respective step folder
    - e.g. `Analysis_Data/<step_folder>/plots/`
    - Create if missing, e.g. `mkdir -p Analysis_Data/<step_folder>/plots/`
    - <step_folder> is the folder where the step is executed (e.g. `alignment`, `counts`, `deeptools`, etc.)
- Use Python with matplotlib/seaborn/pandas — these are standard in bioinformatics envs
- Always save figures as `.svg` and print the output path
- After saving, describe what the figure shows and flag any sample that looks anomalous
- Never display figures interactively (no `plt.show()`) — save only

---

## Check environment first
```bash
python3 -c "import matplotlib, seaborn, pandas; print('OK')" 2>&1
```
If a library is missing:
- Update the `environment.yml` file with the missing libraries
- Update environment: `conda env update --file environment.yml`

---

## Plot Templates

### 1. Alignment Rate Bar Chart (after step_2)

```python
#!/usr/bin/env python3
"""Plot HISAT2 alignment rates from step_2 .error logs."""
import os, re, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

LOG_DIR = "Logs"
FIG_DIR = "*/*/plots"
os.makedirs(FIG_DIR, exist_ok=True)

# Parse alignment rates from HISAT2 stderr (captured in .error files)
pattern = re.compile(r"(\S+).*?(\d+\.\d+)% overall alignment rate", re.S)
records = []
for fname in sorted(os.listdir(LOG_DIR)):
    if fname.endswith(".error") and "step_2" in fname:
        text = open(os.path.join(LOG_DIR, fname)).read()
        m = re.search(r"(\d+\.\d+)% overall alignment rate", text)
        if m:
            sample = fname.replace(".error", "")
            records.append({"sample": sample, "alignment_rate": float(m.group(1))})

if not records:
    print("No alignment rates found in log/. Check step_2 .error files.")
    exit(1)

df = pd.DataFrame(records).sort_values("alignment_rate")

fig, ax = plt.subplots(figsize=(max(8, len(df)*0.5), 5))
colors = ["#d9534f" if r < 70 else "#5cb85c" for r in df["alignment_rate"]]
ax.barh(df["sample"], df["alignment_rate"], color=colors)
ax.axvline(70, color="orange", linestyle="--", label="70% threshold")
ax.axvline(50, color="red", linestyle="--", label="50% (HALT)")
ax.set_xlabel("Overall Alignment Rate (%)")
ax.set_title("HISAT2 Alignment Rates — <PROJECT_NAME>")
ax.legend()
plt.tight_layout()

out = os.path.join(FIG_DIR, "step2_alignment_rates.svg")
plt.savefig(out, dpi=300)
print(f"Saved: {out}")

# Flag problem samples
problems = df[df["alignment_rate"] < 70]
if not problems.empty:
    print("\n⚠️  SAMPLES BELOW 70% ALIGNMENT RATE:")
    print(problems.to_string(index=False))
```

---

### 2. Read Count Distribution (after step_4)

```python
#!/usr/bin/env python3
"""Plot read count distributions per sample from featureCounts output."""
import os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import glob

FIG_DIR = "plots"
os.makedirs(FIG_DIR, exist_ok=True)

# Adjust glob to match your featureCounts output file
count_files = glob.glob("Counts/*.counts") or glob.glob("Counts/*.txt")
if not count_files:
    print("No count files found in Counts/. Adjust path.")
    exit(1)

dfs = []
for f in sorted(count_files):
    sample = os.path.basename(f).replace(".counts","").replace(".txt","")
    df = pd.read_csv(f, sep="\t", comment="#", index_col=0)
    count_col = df.columns[-1]   # featureCounts: last column is counts
    dfs.append(df[count_col].rename(sample))

counts = pd.concat(dfs, axis=1).fillna(0)
log_counts = counts.apply(lambda x: (x+1).apply(__import__('math').log2))

fig, ax = plt.subplots(figsize=(max(8, len(count_files)*0.6), 5))
log_counts.boxplot(ax=ax, rot=45)
ax.set_ylabel("log2(counts + 1)")
ax.set_title("Read Count Distribution per Sample")
plt.tight_layout()

out = os.path.join(FIG_DIR, "step4_count_distribution.svg")
plt.savefig(out, dpi=300)
print(f"Saved: {out}")

# Flag samples with very low total counts
totals = counts.sum().sort_values()
low = totals[totals < totals.median() * 0.3]
if not low.empty:
    print("\n⚠️  SAMPLES WITH VERY LOW TOTAL COUNTS:")
    print(low.to_string())
```

---

### 3. FastQC Summary Table (after step_0f)

```python
#!/usr/bin/env python3
"""Summarize FastQC pass/warn/fail across samples."""
import os, glob, zipfile, re, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

FIG_DIR = "plots"
os.makedirs(FIG_DIR, exist_ok=True)

qc_zips = glob.glob("FastQC/**/*.zip", recursive=True)
records = []
for zpath in qc_zips:
    sample = os.path.basename(zpath).replace("_fastqc.zip","")
    with zipfile.ZipFile(zpath) as z:
        summary_name = [n for n in z.namelist() if n.endswith("summary.txt")][0]
        with z.open(summary_name) as f:
            for line in f.read().decode().splitlines():
                status, module, _ = line.split("\t")
                records.append({"sample": sample, "module": module, "status": status})

if not records:
    print("No FastQC zip files found. Check FastQC output directory.")
    exit(1)

df = pd.DataFrame(records)
pivot = df.pivot_table(index="module", columns="sample", values="status", aggfunc="first")

color_map = {"PASS": "#5cb85c", "WARN": "#f0ad4e", "FAIL": "#d9534f"}
cell_colors = pivot.applymap(lambda x: color_map.get(x, "white"))

fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns)*0.8), max(6, len(pivot)*0.4)))
ax.axis("off")
tbl = ax.table(cellText=pivot.values, rowLabels=pivot.index,
               colLabels=pivot.columns,
               cellColours=cell_colors.values, loc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
ax.set_title("FastQC Summary", fontsize=12, pad=20)
plt.tight_layout()

out = os.path.join(FIG_DIR, "step0f_fastqc_summary.svg")
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
```

---

## Interpretation Guide

After generating a plot, always provide a brief written interpretation:

**Alignment rates:**
- >90% = excellent
- 70–90% = acceptable
- 50–70% = investigate (wrong genome? adapter contamination? wrong strand?)
- <50% = HALT — report to user before proceeding

**Count distributions:**
- Boxes should be roughly at the same median level across samples
- Outlier samples (median very low or very high) may need to be flagged or excluded
- Presence of many zeros may indicate alignment or counting issues

**FastQC:**
- `Per base sequence quality` FAIL → trimming may have failed
- `Adapter Content` FAIL → trimming may have failed (check step_1)
- `Per sequence GC content` FAIL → possible contamination or wrong species
- `Sequence Duplication Levels` WARN/FAIL → expected in RNA-seq and enrichment-based protocols (ChIP/CUT&RUN), usually acceptable
