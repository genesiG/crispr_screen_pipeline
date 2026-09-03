# SKILL: skill_log_inspector.md
# PURPOSE: Read and classify LSF .error/.log files after any pipeline job. Routes to skill_troubleshoot.md on failure.

## Step 1 — Find logs

Log files land in `log/` (or subfolders) of the step's working dir:
```bash
find . -name "*.error" -newer Scripts/step_N_<name>.py | sort
ls -lt log/*.error 2>/dev/null | head -20
```

## Step 2 — Triage (run in order; stop at first hit)

```bash
LOG="log/<sample>.error"
tail -50 "$LOG"                         # always start here; never cat without wc -l first
grep -n -iE "traceback|TERM_MEMLIMIT|oom-kill|killed|filenotfound|no such file|cannot open" "$LOG" | head -20
grep -n -iE "error|fail|abort|exception" "$LOG" | grep -v "^#" | head -20
grep -n -iE "warn|skip|not found|missing" "$LOG" | head -10
```
> Never `cat` logs >1000 lines. Run `wc -l "$LOG"` first; use `head -20` + `tail -80` + keyword grep.

## Step 3 — Classify → action

| Class | Signature | Action |
|---|---|---|
| ✅ CLEAN | No error keywords; outputs non-empty | Advance to next step |
| ⚠️ WARN | `warn`/`skip` only; outputs exist | Note warnings, advance |
| ❌ PYTHON | `Traceback (most recent call last)` | Load **skill_troubleshoot.md §A** |
| ❌ TOOL CRASH | Tool-native error message, no traceback | Load **skill_troubleshoot.md §B** |
| ❌ OOM/TIMEOUT | `TERM_MEMLIMIT` / `oom-kill` / `Killed` | Load **skill_troubleshoot.md §C** |
| ❌ PATH | `FileNotFoundError` / `No such file` / `cannot open` | Load **skill_troubleshoot.md §D** |
| ❓ UNCLEAR | None of the above match | Load **skill_troubleshoot.md §E** |

## Step 4 — Verify outputs (always, even for CLEAN logs)

Empty outputs = failure, regardless of log status. Check the file type(s) expected from this step:
```bash
# Adjust glob to match step output type:
# BAM: *.bam | trimmed reads: *.fastq.gz | peaks: *.narrowPeak *.broadPeak *.bed
# counts: *.counts *.tsv | bigWig: *.bw | R objects: *.rds *.RData
find . -name "*.bam" -size 0                                    # zero-size = bad
find . -name "*.bam" -newer Scripts/step_N_<name>.py -exec ls -lh {} \;
```

## Step 5 — Alignment rate check (alignment steps only)

```bash
# Covers bowtie2, HISAT2, STAR, bowtie
grep -E "overall alignment rate|Uniquely mapped reads %|at least one reported alignment|aligned concordantly" log/*.error | sort
```
- < 70% → flag for inspection
- < 50% → **HALT, report to user**
