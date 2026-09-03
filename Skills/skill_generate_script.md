# SKILL: skill_generate_script.md
# PURPOSE: Create or update pipeline step scripts matching existing conventions.

## Step 1 — Read source files first (required)
```bash
view Scripts/config.py               # definition of experiment type and other contants
view Skills/skill_reference.md       # understand reference tools for each experiment and task
view Scripts/step_3_process_bam.py   # canonical template — mirror its structure exactly
view Scripts/ldsample.py             # load_samples(), SAMPLES, SAMPLES_CTL
```
Also check `Metadata/sample_metadata.csv` for dataset-specific edge cases.

## Step 2 — Script structure (follow step_3 exactly)

Section order: shebang → docstring → imports → module-level dirs/constants → process_function() → main() → `if __name__ == "__main__"`.

**Key constraints** (deviating from any of these is a bug):

| What | Rule |
|---|---|
| Paths | 100% from `config.*` — no hardcoded strings |
| `log_dir` | Always `os.path.join(work_dir, "log")` |
| Batch + log files | Always written inside `log_dir/` |
| Job name | `f"{sample_name}.<stepverb>"` (no spaces) |
| `conda activate` | Inside the batch heredoc, never before `bsub` |
| `#BSUB -P` | Always `config.PROJECT_NAME` |
| Resource footer | `bjobs -l $LSB_JOBID | grep -A10 "Resource usage"` — last line of every batch |
| Submission | `os.system(f"bsub < {batch_file}")` then `time.sleep(1)` |
| PE/SE branch | `main()` must mirror step_3's `IS_PAIRED_END` pattern exactly |
| Guard | `if __name__ == "__main__": main()` always present |

**Conda env:** check the project's `environment.yml`.

## Step 3 — config.py addition (required for every new step)

Append to the relevant section:
```python
# Step N — <description>
STEP_N_OUT_DIR = os.path.join(<PARENT_DIR>, "<subfolder>")
```
Then `os.makedirs(config.STEP_N_OUT_DIR, exist_ok=True)` in the script.

## Step 4 — Validate before delivering
```bash
python3 -m py_compile Scripts/step_N_*.py && echo "OK"
grep -c "time.sleep\|bsub\|log_dir\|conda activate" Scripts/step_N_*.py
```

## NOTE on R — reticulate integration

When building R scripts, always call `use_python(Sys.which("python"), required = TRUE)` **before** any `py_run_file("config.py")` or `py$*` access. Failure to do so causes silent fallback to system Python and wrong path resolution.

```r
# Example for ChIP-seq, CUT and RUN, and CUT and TAG pipelines
library(reticulate, quietly = TRUE, verbose = FALSE)
use_python(Sys.which("python"), required = TRUE)
py_run_file("Scripts/config.py")
# Now py$WORKDIR, py$BAMDIR etc. are available
```

---