# Project Rules & Automatic GitHub Synchronization

This repository is configured to synchronize automatically with GitHub (`https://github.com/vmcerotres/report-2-dimension-reduction-methods.git`) whenever files are created or modified.

## Mandatory Agent Workflow

1. **Before making edits in a new session:** Ensure the latest changes from collaborators are pulled (`git pull --rebase --autostash origin main`).
2. **LaTeX Source in Root, PDF Compiled into `report_2/`:**
   - Edit `Report_2_PCA.tex` in the workspace root whenever report changes are requested.
   - Always compile `Report_2_PCA.tex` to `report_2/Report_2_PCA.pdf` (removing `.aux`, `.log`, `.out` afterwards) so `report_2/Report_2_PCA.pdf` is always up to date.
3. **Clean Deliverable Directory (`report_2/`):**
   - `report_2/` contains only `script.py`, `Report_2_PCA.pdf`, the 8 PNG plots, and `forty+soybean+cultivars+from+subsequent+harvests/`.
4. **Every time any file is created, modified, or regenerated:**
   - The lifecycle hooks in `.agents/hooks.json` automatically run `.agents/auto_sync.sh` to recompile the PDF (if needed) and commit/push all changes to GitHub.
