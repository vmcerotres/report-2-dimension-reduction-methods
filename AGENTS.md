# Project Rules & Automatic GitHub Synchronization

This repository is configured to synchronize automatically with GitHub (`https://github.com/vmcerotres/report-2-dimension-reduction-methods.git`) whenever files are created or modified.

## Mandatory Agent Workflow

1. **Before making edits in a new session:** Run `git pull --rebase --autostash origin main` (or allow the `SessionStart` hook to pull the latest changes from collaborators).
2. **Every time any file is created, modified, or regenerated:**
   - The `.agents/hooks.json` (`PostToolUse` and `Stop` hooks) automatically executes `.agents/auto_sync.sh` to stage, commit, and push all changes to `origin main`.
   - If you run a script or modify files directly and want to ensure immediate synchronization, run:
     ```bash
     ./.agents/auto_sync.sh
     ```
3. **Keep Plots, Report, and `report_2/` Synchronized:**
   - The clean deliverable folder for the professor is `report_2/` (which must remain completely free of any `.agent`, `.agents`, `.git`, or AI-related references).
   - Whenever `script.py` or `Report_2_PCA.md` is modified, always keep `report_2/script.py`, `report_2/Report_2_PCA.md`, and the PNG charts inside `report_2/` updated and synchronized.
