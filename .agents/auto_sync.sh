#!/bin/sh
# Auto-sync script executed by Antigravity hooks (.agents/hooks.json) and/or fswatch watcher.
# Always outputs valid empty JSON `{}` on stdout so it satisfies PostToolUse/SessionStart hook contracts.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT" || { echo "{}"; exit 0; }

# Consume stdin if invoked as an Antigravity hook
if [ ! -t 0 ]; then
  cat >/dev/null 2>&1
fi

{
  if [ -d ".git" ]; then
    # Pull any remote updates quietly
    git pull --rebase --autostash origin main >/dev/null 2>&1 || true

    # Check if there are local changes (modified, added, deleted, or untracked)
    if [ -n "$(git status --porcelain 2>/dev/null)" ]; then
      git add -A >/dev/null 2>&1
      CHANGED_FILES="$(git diff --cached --name-only | tr '\n' ' ' | sed 's/ *$//')"
      TIMESTAMP="$(date '+%Y-%m-%d %H:%M:%S')"
      git commit -m "Auto-sync ($TIMESTAMP): $CHANGED_FILES" >/dev/null 2>&1 || true
      git push -u origin main >/dev/null 2>&1 || true
    fi
  fi
} >/dev/null 2>&1

echo "{}"
