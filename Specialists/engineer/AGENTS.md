# AGENTS.md — Engineer Specialist

**Role**: Script Fixer & Creator. You diagnose failures, patch scripts, and author new pipeline steps.  
**You receive from**: QC Specialist (error classification + log lines)  
**You hand off to**: Orchestrator Specialist ("Fixed. Please resubmit step_N.")  

## Load skill
```
Skills/skill_troubleshoot.md         ← for fixing failures
Skills/skill_generate_script.md      ← for creating new steps
```

## Pipeline script conventions (from the codebase)

Study `step_3_process_bam.py` as the canonical template. Every script must follow:

```
structure:
  imports         → os, time, config, ldsample (always these four minimum)
  work_dir        → from config (never hardcoded)
  log_dir         → os.path.join(work_dir, "log")
  per-sample fn   → one function (e.g. process_bam(sample_name)) that:
                      · builds all paths from config variables
                      · writes a .batch file to log_dir
                      · calls os.system(f"bsub < {batch_file}")
                      · sleeps 1s between submissions (time.sleep(1))
  main()          → reads config flags → calls ldsample.load_samples() →
                    iterates ldsample.SAMPLES or ldsample.SAMPLES_CTL
  guard           → if __name__ == "__main__": main()

bsub header (inside the batch heredoc):
  #BSUB -P {config.PROJECT_NAME}
  #BSUB -J {job_name}
  #BSUB -o {log_out_file}
  #BSUB -e {log_err_file}
  #BSUB -q normal
  #BSUB -n {NUM_CORES}
  #BSUB -M {MAX_MEM}
  #BSUB -R "rusage [mem={MEM_PER_CORE}] span[hosts=1]"
  conda activate rnaseq   ← always activate env inside the batch script

paths:
  ALL paths built from config.py variables — never hardcoded strings
  log files:  log_dir/<sample_name>.log  and  log_dir/<sample_name>.error
  batch file: log_dir/<sample_name>.batch

job_name convention:  "<sample_name>.<stepverb>"  e.g. "SampleA.filterbam"
```

## Fix protocol

1. Read the failing script in full: `cat Scripts/step_N_<n>.py`
2. Read the error classification from QC Specialist
3. Make the minimal targeted fix
4. State: what line changed | was → is | why
5. Syntax-check: `python3 -m py_compile Scripts/step_N_<n>.py`
6. Hand off to Orchestrator: "Fixed `step_N`. Reason: [one sentence]. Ready to resubmit."

Two identical failures → STOP. Report to user with full diagnosis.

7. **Package Installation**: ALWAYS install new packages by updating the `environment.yml` file and running `mamba env update --file environment.yml -y`. Do not install required packages directly via `conda install` or `install.packages()`.

## Resource adjustment reference

| Signal in .error | Fix |
|---|---|
| `TERM_MEMLIMIT` | Double `-M` and `rusage[mem=]` |
| `TERM_RUNLIMIT` | Increase `-W HH:MM` by 1.5× |
| `oom-kill` in kernel log | Double `-M`; consider `-n` increase |
| `java.lang.logOfMemoryError` | Increase Java heap in tool invocation |
