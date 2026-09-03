---
name: skill_tool_registry
description: >
  Reference registry of approved tools, packages, and parameters for every
  major analysis type in the chromatin/transcriptomics pipeline. Load this
  skill whenever you need to select or justify a tool for peak calling
  (ChIP-seq, CUT&RUN, CUT&TAG, ATAC-seq), differential binding, differential
  gene expression, normalization, QC, visualization, or CRISPR screens.
  Also load this skill whenever skill_generate_script.md needs to confirm the
  correct tool stack before writing analysis code. Trigger on any mention of
  "which tool", "peak calling", "differential binding", "DESeq2 vs edgeR",
  "normalization method", "CRISPR screen analysis", "what package should I
  use", or "recommended pipeline" for genomics data.
---

# Tool & Technology Registry — Chromatin Profiling & Transcriptomics Pipeline

This is the authoritative reference for tool selection across all assay types
handled by this pipeline. Downstream skills (especially `skill_generate_script`)
**must** consult this registry before emitting code or bsub scripts.

---

## 0 — How to Use This Skill

| Caller | Action |
|--------|--------|
| `skill_generate_script` | Look up the assay row → read the full section for that analysis type → apply tool + parameter defaults |
| Direct user query ("which peak caller should I use?") | Identify assay → return the recommendation table + rationale blurb |
| Ambiguous assay | Ask one clarifying question: **"Is this ChIP-seq, CUT&RUN/CUT&TAG, or ATAC-seq?"** before proceeding |

---

## 1 — Quick-Reference Master Table

| Analysis | Primary Tool | Language / Env | Key Flag / Note |
|----------|-------------|----------------|-----------------|
| QC — raw reads | FastQC + MultiQC | CLI | Always run pre- and post-trim |
| Adapter trimming | Trimmomatic | CLI | PE: `ILLUMINACLIP`, SE: same |
| Alignment — ChIP/CUT/ATAC | Bowtie2 | CLI | `--very-sensitive --no-mixed --no-discordant` (PE) |
| Alignment — RNA-seq | STAR | CLI | `--outSAMtype BAM SortedByCoordinate` |
| BAM processing | samtools + bedtools | CLI | See §3 |
| Duplicate marking/removal | Picard `MarkDuplicates` | CLI | Pin `picard=2.27.5` in env |
| Peak calling — ChIP-seq | MACS3 | Python/CLI | `--broad` for broad marks; `--nomodel` optional |
| Peak calling — CUT&RUN/CUT&TAG | MACS3 (narrow marks) / SEACR (broad / low-signal) | Python/CLI | See §4 |
| Peak calling — ATAC-seq | MACS3 `--nomodel --shift -100 --extsize 200` | Python/CLI | Fragment-centered |
| Differential binding | csaw (R) | R / Bioconductor | Window-based; handles low-signal assays well |
| Differential gene expression | DESeq2 (default) / edgeR | R / Bioconductor | DESeq2 for ≥3 reps; edgeR for 2 reps |
| Normalization — spike-in | Custom `step_4c_spikein_normalization.py` | Python | EpiCypher CUTANA protocol |
| Normalization — reference-free | ChIPseqSpikeInFree (R) | R | Bin-based; requires ≥4 samples |
| Normalization — count-based | TMM via edgeR (`calcNormFactors`) | R | For RNA-seq or ChIP without spike-in |
| Coverage tracks | deepTools `bamCoverage` / `bamCompare` | Python/CLI | RPGC or custom SF from norm step |
| Heatmaps / profiles | deepTools `computeMatrix` + `plotHeatmap` / EnrichedHeatmap (R) | CLI / R | See §8 |
| Clustering | K-means (Python `sklearn`) + silhouette selection | Python | See §9 |
| CRISPR screen | MAGeCK (`mageck count` + `mageck test`) | Python/CLI | MLE for complex designs |
| Motif analysis | HOMER `findMotifsGenome.pl` | CLI | Needs genome FASTA |
| Gene annotation of peaks | ChIPseeker (R) | R / Bioconductor | `annotatePeak()` with Bioconductor TxDb |

---

## 2 — Environment & Infrastructure

```
Conda environment : chipseq (or rnaseq or crisprscreen)
Scheduler         : LSF / bsub
Pipeline language : Python (job submission) + R (statistical analysis)
Config module     : config.py  (WORKDIR, BAMDIR, PEAKDIR, SCALINGDIR, etc.)
Sample metadata   : ldsample.py / samples.txt / sample_metadata.txt
```

## 3 — BAM Processing Standards

| Step | Tool | Parameters / Notes |
|------|------|--------------------|
| Sort | `samtools sort` | `-@ 4` (threads plateau here) |
| Index | `samtools index` | Required after every sort or filter |
| Mark duplicates | `picard MarkDuplicates` | `picard=2.27.5` pinned in env; `REMOVE_DUPLICATES=true` when `config.REMOVE_DUPLICATES=True` |
| MAPQ filter | `samtools view -q {MAPQ}` | Default `MAPQ=10` from `config.py` |
| Blacklist removal | `bedtools intersect -v -abam` | Use `config.BLACKLIST`; do **not** convert to BED — preserve SAM tags |
| Spike-in split | `samtools view` on contig names | Must occur **before** human-genome filters; see §6 |

**Anti-pattern:** BAM → BED → BAM round-trips lose SAM auxiliary tags and FLAG values. Always stay BAM-native with `samtools` + `bedtools intersect`.

---

## 4 — Peak Calling

### 4a ChIP-seq

```
Tool      : MACS3
Genome    : config.GENOME_SIZE_FOR_MACS  (3.1e9 for t2t, 2.7e9 for mm39, "hs" otherwise)
Broad marks (H3K27me3, H3K9me3, H3K27ac domain): --broad --broad-cutoff 0.1 -q 0.01
Narrow marks (H3K4me3, CTCF, TF): -q 0.01  (no --broad)
Input control : supply with -c when available; omit for IgG-only designs
Key config    : config.MACS_Q_VALUE, config.MACS_CALL_BROAD_PEAKS, config.MACS_MAX_GAP
```

### 4b CUT&RUN / CUT&TAG

CUT&RUN/CUT&TAG produces shorter, sharper fragments than ChIP-seq. Tool choice depends on signal character:

| Mark type | Signal level | Recommended caller | Mode |
|-----------|-------------|--------------------|------|
| Narrow / transcription factor | High | MACS3 | `--nomodel -q 0.01` |
| Broad histone (H3K27me3, H3K9me3) | High | MACS3 | `--broad --nomodel` |
| Any mark | Low / noisy | SEACR | `non stringent` mode, IgG as control |

SEACR is preferred when signal-to-noise is poor or sample depth is low (<5 M mapped reads after dedup). MACS3 is preferred otherwise for consistency with the rest of the pipeline.

**Fragment size filter:** For CUT&RUN, pre-filter BAM to ≤200 bp fragments (nucleosome-free) or 180–700 bp (mono-nucleosome) depending on target, using `samtools view -F 4` and awk on TLEN.

### 4c ATAC-seq

```
Tool    : MACS3
Flags   : --nomodel --shift -100 --extsize 200 --nolambda -q 0.05
Rationale: Centers signal on the Tn5 cut site; --nolambda avoids background
           inflation in open-chromatin deserts
```

---

## 5 — Differential Binding

```
Tool        : csaw (R / Bioconductor)
Why csaw    : Window-based counting avoids pre-specifying peak boundaries;
              handles broad marks and low-signal CUT&RUN better than
              count-in-peaks approaches (DiffBind)
readParam   : minq = config.MAPQ, pe = "both" (PE) or "first" (SE),
              max.frag = 500, dedup = config.REMOVE_DUPLICATES
Window size : 150 bp (narrow marks / TF), 500–1000 bp (broad histone marks)
Merging     : mergeResults() + combineTests() with FDR ≤ 0.05
Normalization: feed pre-computed spike-in / TMM / ChIPseqSpikeInFree SFs
              into csaw via normOffsets() or directly via sizeFactors()
```

**Parameters reference file:** `csaw_parameters.R` in project Scripts.

---

## 6 — Differential Gene Expression

| Condition | Tool | Rationale |
|-----------|------|-----------|
| ≥3 biological replicates per group | **DESeq2** | Empirical Bayes shrinkage; robust size-factor estimation |
| 2 replicates per group | **edgeR** | GLM handles low-replicate designs better |
| Time-course / multi-factor | DESeq2 with LRT | `test="LRT"`, reduced model drops the factor of interest |
| Single-cell pseudobulk | DESeq2 on aggregated counts | Aggregate per sample, not per cell |

### Multiple-testing correction

Always use `padj` (Benjamini-Hochberg FDR) from DESeq2 / edgeR `p.adjust(method="BH")`. Report both raw p-value and adjusted p-value. Threshold: FDR ≤ 0.05, |log2FC| ≥ 1 (adjust per experiment).

---

## 7 — Visualization (Coverage Tracks & Heatmaps)

### BigWig generation

```
Tool     : deepTools bamCoverage (via step_5_generate_bigwig.py)
Bin size : config.BIN_SIZE (default 10 bp)
Norm     : --scaleFactor {SF} computed in §5; or built-in RPGC for
           exploratory tracks only
```

### Heatmaps / profile plots

| Use case | Tool | Notes |
|----------|------|-------|
| Standard profile / heatmap | deepTools `computeMatrix` + `plotHeatmap` | `reference-point` or `scale-regions` mode |
| Clustering-aware heatmap | EnrichedHeatmap (R) via `step_6d_build_matrix.R` | Requires `normalizeToMatrix(mean_mode="w0")` |
| MA plots | `smoothScatter` + `abline(0,1)` (R) | Used in normalization comparison (`step_4d_compare_normalizations.R`) |
| Complex multi-sample heatmap | ComplexHeatmap (R) | Combine with K-means cluster labels |

**Centering rule:** Always resize peak/window GRanges to 1-bp centers before `normalizeToMatrix()` — `resize(targets, width=1L, fix="center")` — to prevent start-position skew on wide csaw windows.

---

## 8 — Clustering

```
Method   : K-means (scikit-learn KMeans, n_init=50, algorithm="lloyd")
k range  : Sweep k_min..k_max; select best k by silhouette score
Parallelism: joblib Parallel across k values (n_jobs = NUM_CORES)
Input    : Combined signal matrix (samples × bins), NaN → 0
Output   : cluster_assignments_k{k}.csv (1-indexed), silhouette_scores.csv,
           best_k.txt
```

See `step_6d_clustering.py` for the full submit/run pattern.

---

## 9 — CRISPR Screen Analysis

| Step | Tool | Command / Notes |
|------|------|-----------------|
| Count sgRNA reads | `mageck count` | `--list-seq sgrna_library.txt --fastq *.fastq.gz` |
| Differential essentiality | `mageck test` | `--treatment-id` / `--control-id`; outputs gene-level RRA scores |
| Complex designs (multiple conditions, covariates) | `mageck mle` | Requires design matrix; reports beta scores per condition |
| QC — read mapping rate | `mageck count` log | Expect ≥60% mapping; flag samples <50% |
| Visualization | `mageck pathway` or custom R (ggplot2 + MAGeCKFlute) | MAGeCKFlute wraps standard downstream plots |

**Statistical note:** MAGeCK RRA aggregates sgRNA-level p-values by gene using a robust rank aggregation algorithm. Report FDR (`fdr` column) not raw p-value for gene-level hits.

---

## 10 — Anti-Patterns (Never Do These)

| Anti-pattern | Correct approach |
|---|---|
| BAM → BED → BAM for filtering | Use `bedtools intersect -abam` to stay BAM-native |
| Spike-in split after human filters | Split first, filter human separately |
| `use_python()` instead of `use_condaenv()` in R | `use_condaenv("chipseq", required=TRUE)` always |
| Picard without version pin | Pin `picard=2.27.5` in `environment.yml` |
| Bowtie2 for RNA-seq | Use STAR |
| DiffBind for broad histone marks | Use csaw with appropriate window size |
| `plotHeatmap` without `computeMatrix` | Always run `computeMatrix` first; do not reuse stale `.gz` matrices |
| Over-provisioning cores for bedtools/Picard | Both are single-threaded; `samtools` caps at ~4 threads |
| edgeR for ≥3 replicates without justification | Default to DESeq2; document if edgeR is chosen |
| MAGeCK RRA p-values without FDR | Always report `fdr` column as the significance metric |

---

## 11 — Reference Files (load as needed)

| File | Load when |
|------|-----------|
| `environment.yml` | Learn which environment configuration is being used |
| `Scripts/utils/csaw_parameters.R` | Writing any differential binding script |
| `Scripts/step_4c_spikein_normalization.py` | Implementing or debugging spike-in SF |
| `Scripts/step_4d_compare_normalizations.R` | Comparing normalization methods / MA plots |
| `Scripts/step_5_generate_bigwig.py` | Writing bigWig generation jobs |
| `Scripts/step_6d_clustering.py` | Implementing K-means clustering step |
| `Scripts/step_6d_build_matrix.R` | Building `normalizedMatrix` for EnrichedHeatmap |
