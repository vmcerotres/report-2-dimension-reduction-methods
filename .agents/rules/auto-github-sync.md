# Project Rules & Automatic GitHub Synchronization

This repository is configured to synchronize automatically with GitHub (`https://github.com/vmcerotres/report-2-dimension-reduction-methods.git`) whenever files are created or modified.

## Mandatory Agent Workflow

1. **Before making edits in a new session:** Ensure the latest changes from collaborators are pulled (`git pull --rebase --autostash origin main`).
2. **Every time any file is created, modified, or regenerated:**
   - The lifecycle hooks in `.agents/hooks.json` automatically run `.agents/auto_sync.sh` to commit and push all changes to GitHub.
   - If needed, you can also trigger synchronization manually via `./.agents/auto_sync.sh`.
