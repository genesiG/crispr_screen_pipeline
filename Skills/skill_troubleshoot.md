# SKILL: skill_troubleshoot.md
# PURPOSE: Fix diagnosed pipeline failures across any bioinformatics pipeline.
# If called directly without prior log inspection, run skill_log_inspector.md Steps 2–3 first, then jump to the relevant section below.

## Rules (non-negotiable)
1. `view` the full script before any edit — never edit blind.
2. One fix at a time — minimum change only.
3. Explain: what line, what it was, what it is now, why.
4. Resubmit after fix; go back to log inspection.
5. Same error recurs twice → **STOP, report to user**.
6. New packages → update `environment.yml`, run `mamba env update --file environment.yml -y`. Never `conda install` or `install.packages()` directly.

---

## §A — Python Exception

```bash
grep -A 20 "Traceback" log/<sample>.error
view Scripts/<step>.py                      # read before editing
```

| Error | Likely cause | Fix |
|---|---|---|
| `KeyError` on sample dict | Column name mismatch in sample sheet | Check sample loader (`ldsample.py` or equivalent) vs step expectation |
| `FileNotFoundError` | Config path wrong or previous step incomplete | Verify config paths; check upstream outputs |
| `CalledProcessError` | External tool failed inside Python | Read tool stderr lines above the exception |
| `IndexError`/`StopIteration` | Samples file empty or malformed | `head` the samples file; verify delimiter and column names |
| `ModuleNotFoundError` | Package missing from environment | Add to `environment.yml`; run `mamba env update` (Rule 6) |

---

## §B — Tool Crash

```bash
grep -iE "error|exception|fail" log/<sample>.error | grep -v "^#" | head -30
```

| Tool | Common error | Fix |
|---|---|---|
| **Any aligner** (bowtie2/STAR/HISAT2) | `Cannot open index` / index not found | Check index path in config; verify index was built |
| **Any aligner** | Low alignment rate | Check FASTQ quality; verify genome build matches library prep |
| **samtools** | `bgzf_read: Read block failed` | BAM truncated — rerun upstream step |
| **samtools/sambamba/picard** | `Cannot allocate memory` | Increase `-M` in bsub; ensure `span[hosts=1]` |
| **MACS3/SEACR/csaw** | No peaks / empty output | Verify BAM has reads; check input/control pairing |
| **featureCounts/HTSeq** | No features counted / GTF mismatch | Check chr naming between BAM and GTF (`chr1` vs `1`) |
| **MAGeCK** | Library file not found | Check sgRNA library path in config |
| **bedtools** | Chromosome name mismatch | Verify genome sizes file matches BAM contig names |
| **Rscript** | `there is no package called X` | Add to `environment.yml`; run `mamba env update` (Rule 6) |
| **Java tools** (Trimmomatic/Picard) | `OutOfMemoryError` / GC overhead | Add `-Xmx<N>g` to the tool invocation in the batch script |

---

## §C — LSF Resource Kill (`TERM_MEMLIMIT` / OOM / Timeout)

```bash
grep -n "BSUB\|-M \|-W \|-n " Scripts/<step>.py
```

Adjust in the script's `batch_cmd` f-string:
- OOM → double `MAX_MEM`; recalculate `MEM_PER_CORE = MAX_MEM / NUM_CORES`
- Timeout → increase `-W HH:MM` by 1.5×
- Both `-M` and `-R "rusage[mem=...]"` must be updated together.

---

## §D — Path / Config Error

If the pipeline uses a `config.py`:
```bash
python3 -c "
import sys; sys.path.insert(0,'Scripts'); import config, os
for k,v in vars(config).items():
    if not k.startswith('_') and isinstance(v, str) and '/' in v:
        print(k, '=', v, '→', 'OK' if os.path.exists(v) else '*** MISSING ***')
"
```

If no `config.py`, inspect the script directly for path variables:
```bash
grep -n "os.path\|join\|\.bam\|\.fastq\|\.gtf\|\.bed" Scripts/<step>.py | head -30
```

If a path is missing: (1) created by a previous step that didn't finish? (2) typo in config? (3) filesystem remount?

---

## §E — Unknown / Unclear

```bash
wc -l log/<sample>.error
head -30 log/<sample>.error
tail -80 log/<sample>.error
find . -newer Scripts/step_N_<name>.py -type f | sort
```

If still unclear: **STOP** — report to user with last 50 lines of `.error`, the bsub command from the script, and your best hypothesis.
