#!/bin/sh
# Auto-sync script executed by hooks (.agents/hooks.json).
# Automatically compiles Report_2_PCA.tex -> report_2/Report_2_PCA.pdf when needed,
# cleans LaTeX auxiliary files, and commits/pushes all changes to GitHub.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT" || { echo "{}"; exit 0; }

if [ ! -t 0 ]; then
  cat >/dev/null 2>&1
fi

export PATH="/Library/TeX/texbin:/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:$PATH"

{
  # Recompile Report_2_PCA.tex -> report_2/Report_2_PCA.pdf if .tex exists and is newer than PDF (or PDF missing)
  if [ -f "Report_2_PCA.tex" ] && command -v pdflatex >/dev/null 2>&1; then
    if [ ! -f "report_2/Report_2_PCA.pdf" ] || [ "Report_2_PCA.tex" -nt "report_2/Report_2_PCA.pdf" ]; then
      pdflatex -interaction=nonstopmode -halt-on-error -output-directory=report_2 Report_2_PCA.tex >/dev/null 2>&1 || true
      pdflatex -interaction=nonstopmode -halt-on-error -output-directory=report_2 Report_2_PCA.tex >/dev/null 2>&1 || true
      rm -f report_2/*.aux report_2/*.log report_2/*.out report_2/*.toc report_2/*.fls report_2/*.fdb_latexmk report_2/*.synctex.gz
    fi
  fi

  if [ -d ".git" ]; then
    git pull --rebase --autostash origin main >/dev/null 2>&1 || true

    if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
      git add -A >/dev/null 2>&1
      CHANGED_FILES="$(git diff --cached --name-only | tr '\n' ' ' | sed 's/ *$//')"
      TIMESTAMP="$(date '+%Y-%m-%d %H:%M:%S')"
      git commit -m "Update ($TIMESTAMP): $CHANGED_FILES" >/dev/null 2>&1 || true
      git push -u origin main >/dev/null 2>&1 || true
    fi
  fi
} >/dev/null 2>&1

echo "{}"
