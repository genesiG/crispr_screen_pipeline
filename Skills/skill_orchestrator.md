# SKILL: skill_orchestrator.md
# Loaded by: Orchestrator Specialist
# PURPOSE: Step-by-step execution template for running pipeline steps on an LSF-based HPC cluster.
#          Experiment-agnostic: works for any sequencing analysis pipeline.
#          Also used by the /run-pipeline workflow for autonomous multi-step execution.

---

## ⚠️ CRITICAL: run_command tool behavior (read before anything else)

Due to tool regressions in the current Antigravity IDE version (e.g. background tasks reporting stale statuses, inconsistent tool-returned log paths, and automatic backgrounding):
**YOU MUST USE A PERSISTENT TERMINAL PATTERN TO AVOID HANGS AND MISSING OUTPUT.**

### The Persistent Terminal & Log Pattern:

1. **First Command in Session**: Always run the first command with `RunPersistent=true` and a large `WaitMsBeforeAsync` (e.g., `30000` or 30 seconds). Keep note of the returned `TerminalID` in the tool output.
2. **Subsequent Commands**: For ALL subsequent `run_command` calls in the session, you **MUST**:
   - Set `RunPersistent` to `true`.
   - Set `RequestedTerminalID` to the `TerminalID` you saved from the first command.
3. **Command Execution & Inline Output**:
   - For **fast commands** (e.g., `hostname`, `whoami`, `ls`, `bjobs`): Set `WaitMsBeforeAsync` to at least `30000` (30 seconds) or more. This ensures the output is returned inline in the tool call response rather than backgrounded.
   - For **long-running commands** (e.g., Python scripts, `bsub` submissions): Set `WaitMsBeforeAsync` to `60000` (60 seconds) or more.
4. **NEVER Use Tool-returned Log Paths, `manage_task`, or `command_status`**:
   - Do NOT call `manage_task status`, `manage_task list`, or `command_status`. They are completely broken/unreliable (e.g. they report `RUNNING` indefinitely).
   - Do NOT try to read the log files in `.system_generated/tasks/` or the `logUri` returned by the tool. They do not reliably correspond to actual command output or files on disk.
5. **Read On-Disk Logs Directly**:
   - Pipeline scripts automatically write their output directly to log files on disk under the `Logs/` directory (e.g., `Logs/step_X_submit.log`, `Logs/step_X_poll.log`, `Logs/step_X_run.log`).
   - Monitor job progress, check submission status, and verify command execution by using the `view_file` tool on the **actual files in the `Logs/` directory** rather than tool-returned task logs.

---

## Pre-session safety check (ALWAYS first)

```bash
hostname   # MUST be nodeXXX.hpc.local
whoami     # must be <your_username>
```

If `hostname` returns anything containing `login`, `mercury`, `consign`, `rhel`, or `hpclogin`:
**Run `bsub -Is bash` from the agent's own terminal** to acquire a compute node before proceeding.
The user's interactive session is a *separate* terminal — the agent cannot reuse it.
Wait for the shell prompt to change to a `nodeXXX` hostname before continuing.

---

## Step Discovery (experiment-agnostic)

Do NOT rely on hardcoded step tables. Discover steps dynamically:

1. **List all `Scripts/step_*.py` files** — the numeric prefix determines execution order
   (e.g., `step_0a` → `step_0b` → `step_1` → `step_2` → `step_2b` → `step_3` → ...)
```bash
ls Scripts/step_*.py | sort
```

2. **Read `config.py`** to determine which steps to skip:
```bash
python3 -c "
import sys; sys.path.insert(0,'Scripts'); import config
print('PROJECT:', config.PROJECT_NAME)
print('EXPERIMENT:', config.EXPERIMENT)
print('USE_TRIMMOMATIC:', config.USE_TRIMMOMATIC)
print('IS_RAW_DATA_FASTQ:', config.IS_RAW_DATA_FASTQ)
print('IS_PAIRED_END:', config.IS_PAIRED_END)
print('REMOVE_DUPLICATES:', config.REMOVE_DUPLICATES)
"
```
   - `USE_TRIMMOMATIC == False` → skip any trimmomatic step
   - `IS_RAW_DATA_FASTQ == True` → skip bam-to-fastq conversion
   - Check `EXPERIMENT` to understand context (rnaseq, chipseq, atacseq, etc.)

3. **Check for `.done` sentinels** — skip steps that already completed:
```bash
ls Logs/*.done 2>/dev/null || echo "No sentinels yet"
```

4. **Determine if step is local or batch**:
   - Read the step script; if it calls `bsub` → batch step (needs polling)
   - If it runs directly → local step (no polling needed)

---

## Executing a Local Step (no bsub)

```bash
python3 Scripts/step_X_<name>.py 2>&1 | tee Logs/step_X_run.log
echo "Exit code: $?"
```
Check exit code and stdout for errors. Hand off to QC Specialist immediately after.

---

## Executing a Batch Step (submits bsub jobs)

### 1. Submit and capture job IDs

```bash
python3 Scripts/step_N_<name>.py 2>&1 | tee Logs/step_N_submit.log

# Capture ALL job IDs — one per sample
mapfile -t JOB_IDS < <(grep -oP '(?<=Job <)\d+(?=>)' Logs/step_N_submit.log)
echo "Submitted ${#JOB_IDS[@]} jobs: ${JOB_IDS[*]}"
```

### 2. Wait — Exponential-Backoff Polling Loop

Run `Scripts/poll_jobs.sh` in the persistent terminal, monitor it by reading `Logs/step_N_poll.log` directly, and wait using the `schedule` tool.

**Algorithm:**
```
timeout = 600          # start at 10 minutes (600 seconds)
MAX_TIMEOUT = 4800     # cap at 80 minutes
MAX_TOTAL = 72000      # hard stop at 20 hours
total_elapsed = 0

LOOP:
  1. Launch the poll script in the persistent terminal:
     run_command("bash Scripts/poll_jobs.sh step_N JOB_IDS --timeout $timeout 2>&1 | tee Logs/step_N_poll.log")
     with RunPersistent=true, RequestedTerminalID=<ID>, WaitMsBeforeAsync=5000.
     This starts the script and backgrounds the IDE task immediately.

  2. Wait for the timeout duration:
     Use the schedule tool to set a timer for the timeout duration (e.g. 600 seconds).
     While waiting, you can periodically inspect Logs/step_N_poll.log using view_file to check progress.

  3. Check results in Logs/step_N_poll.log:
     - If the log ends with "ALL_DONE: X/X jobs completed ✓": Proceed to step 3 (QC).
     - If the log ends with "TIMEOUT: X jobs still running...":
       total_elapsed += timeout
       if total_elapsed >= MAX_TOTAL:
         ABORT: "Jobs exceeded 20h limit"
       timeout = min(timeout * 2, MAX_TIMEOUT)
       Go back to LOOP step 1, but append the log (use `tee -a Logs/step_N_poll.log`).
     - If the log shows failures or exit status 1:
       Proceed to the Error Handling section.
```

**Backoff schedule:**

| Round | Timeout | Cumulative | schedule tool timer | Check log file |
|-------|---------|------------|---------------------|----------------|
| 1     | 10 min  | 10 min     | 10 min (600s)       | Logs/step_N_poll.log |
| 2     | 20 min  | 30 min     | 20 min (1200s)      | Logs/step_N_poll.log |
| 3     | 40 min  | 70 min     | 40 min (2400s)      | Logs/step_N_poll.log |
| 4     | 80 min  | 150 min    | 80 min (4800s)      | Logs/step_N_poll.log |
| 5     | 80 min  | 230 min    | 80 min (4800s)      | Logs/step_N_poll.log |
| ...   | 80 min  | +80 each   | 80 min (4800s)      | Logs/step_N_poll.log |

### 3. QC — Hand off to QC Specialist

After every batch step completes, inspect outputs before advancing:

```
→ QC Specialist: step_N is done.
  Log dir: Logs/  Log prefix: step_N
  Expected outputs: <describe based on step>
  Please inspect and issue your report.
```

Adapt QC checks based on what the step produces (alignment logs, BAMs, count files, peaks, etc.):
```bash
# Example: Check for empty output files
find <OUTPUT_DIR> -name "*.EXT" -size 0
# Example: Check alignment rates from stderr logs
grep "overall alignment rate" Logs/step_N_*.error
```

Act on the QC report before advancing to the next step.

---

## Error Handling

If `poll_jobs.sh` exits 1 (job failure):
1. Check which jobs failed:
```bash
cat Logs/<step_name>.failed
```
2. Inspect the `.error` log for the failed sample(s)
3. If the same step has failed **twice**: STOP and report to the user. Do not loop.
4. Otherwise: diagnose, fix, and retry the step.

---

## After All Steps Complete

Report to the user:
1. Which steps were run vs skipped
2. Any QC warnings (low alignment rates, empty files, etc.)
3. Location of final outputs
