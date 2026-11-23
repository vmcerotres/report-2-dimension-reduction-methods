#!/bin/sh
# Auto-sync script executed by hooks (.agents/hooks.json).
# Also keeps `report_2/` synchronized before committing and pushing to GitHub.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT" || { echo "{}"; exit 0; }

if [ ! -t 0 ]; then
  cat >/dev/null 2>&1
fi

{
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
