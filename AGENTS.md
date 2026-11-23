# Project Rules & Automatic GitHub Synchronization

This repository is configured to synchronize automatically with GitHub (`https://github.com/vmcerotres/report-2-dimension-reduction-methods.git`) whenever files are created or modified.

## Mandatory Agent Workflow

1. **Before making edits in a new session:** Run `git pull --rebase --autostash origin main` (or allow the `SessionStart` hook to pull the latest changes from collaborators).
2. **LaTeX Source Outside, Compiled PDF Inside `report_2/`:**
   - The LaTeX report source file is `Report_2_PCA.tex` in the **workspace root** (outside `report_2/`).
   - Every time any modification is made to the report, data, or code in the chat, `Report_2_PCA.tex` must be compiled into `report_2/Report_2_PCA.pdf` (using `/Library/TeX/texbin/pdflatex`), and all LaTeX auxiliary files (`.aux`, `.log`, `.out`) must be cleaned up so `report_2/` only contains the clean final `.pdf`.
   - The `.agents/auto_sync.sh` hook script also automatically recompiles `Report_2_PCA.tex` -> `report_2/Report_2_PCA.pdf` whenever `Report_2_PCA.tex` or `report_2/script.py` changes.
3. **Single Deliverable Folder (`report_2/`):**
   - All final deliverables for the professor (`script.py`, `Report_2_PCA.pdf`, the 8 PNG charts, and the dataset folder `forty+soybean+cultivars+from+subsequent+harvests/`) live **exclusively** inside `report_2/`. Do not duplicate them in the workspace root, and do not put `.tex` or `.md` files inside `report_2/`.
   - `report_2/` must remain completely clean of any `.agent`, `.agents`, `.git`, auxiliary LaTeX files, or AI-related references so it can be submitted directly to the professor.
4. **Automatic GitHub Push:**
   - The `.agents/hooks.json` (`PostToolUse` and `SessionStart` hooks) automatically executes `.agents/auto_sync.sh` to compile the PDF and stage, commit, and push all changes to `origin main`.
   - If you run a script or modify files directly and want to ensure immediate synchronization, run:
     ```bash
     ./.agents/auto_sync.sh
     ```
