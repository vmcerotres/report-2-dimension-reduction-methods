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
3. **Single Deliverable Folder (`report_2/`):**
   - All project deliverables (`script.py`, `Report_2_PCA.md`, the PNG charts, and the dataset folder `forty+soybean+cultivars+from+subsequent+harvests/`) live **exclusively** inside `report_2/`. Do not duplicate them in the workspace root.
   - `report_2/` must remain completely clean of any `.agent`, `.agents`, `.git`, or AI-related references so it can be submitted directly to the professor.
