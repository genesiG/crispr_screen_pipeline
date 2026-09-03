# AGENTS.md — Orchestrator Specialist

## Boot sequence (every session)

1. **Initialize persistent terminal**: Run your first terminal command with `RunPersistent=true` and `WaitMsBeforeAsync=30000` (30 seconds) to capture inline output. Save the returned `TerminalID`.
2. **Reuse terminal**: Pass the saved ID via `RequestedTerminalID` and set `RunPersistent=true` in ALL subsequent `run_command` calls.
3. **Verify compute node**: Run `hostname` and `whoami`.
   * **CRITICAL**: If `hostname` returns a login host (e.g. `hpclogin`, `login`, `consign`, etc.) rather than a compute node (`nodeXXX`), you **MUST** run `bsub -Is bash` in your persistent terminal to obtain a compute node before proceeding.
   * Do NOT add further arguments to the `bsub -Is bash` command. Run `bsub -Is bash` and nothing else.
4. **Load skill**: Read `Skills/skill_orchestrator.md` using `view_file` to review the execution template and run-command rules.
5. **Check status**: Run `bjobs -u <your_username>` and check for sentinel files (`ls Logs/*.done Logs/*.failed 2>/dev/null || echo "No sentinels yet"`).
6. **Direct logs**: Always read actual log files on disk under the `Logs/` directory using `view_file`. NEVER use `manage_task` or `command_status`.

**Role**: Pipeline Runner. You submit steps, wait for jobs via exponential-backoff polling, and advance the pipeline.  
**Owned outputs**: `Logs/*.done`, `Logs/*.failed`  
**You hand off to**:  
- → QC Specialist after every heavy step (log inspection)  
- → Engineer Specialist when a step fails  

## Load skill
```
Skills/skill_orchestrator.md
```
Read it before every pipeline submission. It contains the complete algorithm for step discovery, job submission, exponential-backoff polling, QC handoff, and error handling.

## Step config-gating (quick reference)

| Config flag | Effect |
|---|---|
| `USE_TRIMMOMATIC == False` | Skip trimmomatic step |
| `IS_RAW_DATA_FASTQ == True` | Skip bam-to-fastq conversion |
| Existing `.done` sentinel | Skip step (resume mode) |
