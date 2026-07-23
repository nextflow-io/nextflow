# Claude Code Hooks for Nextflow Development

This directory contains [Claude Code hooks](https://docs.claude.com/en/docs/claude-code/hooks)
configured to improve the Nextflow development experience. Each hook is a small,
self-contained Bash script under `hooks/`, wired up in `settings.json`.

Hooks are added incrementally — one per pull request — so each can be reviewed and
understood on its own. This file documents the hooks that are currently wired up.

## Requirements

- [`jq`](https://jqlang.github.io/jq/) — hooks parse the Claude Code JSON payload with it.
  If `jq` is not on `PATH`, the hooks no-op silently.

## Hooks

### 1. EditorConfig enforcement (PostToolUse, async)
- **Script**: `hooks/format-editorconfig.sh`
- **Trigger**: after every `Edit`/`Write`/`MultiEdit` on a source file
  (`.groovy`, `.java`, `.gradle`, `.md`, `.txt`, `.yml`, `.yaml`, `.json`)
- **Action**: applies `.editorconfig` rules with `eclint fix`
- **Dependencies**: `eclint` (optional — install with `npm install -g eclint`; the hook
  skips silently if it is not present)
- **`async: true`**: runs in the background so formatting never blocks the conversation

## Enabling / disabling

Hooks are configured in `.claude/settings.json`. Disable one by removing its entry;
adjust timeouts there as needed.

## Troubleshooting

1. Ensure scripts are executable: `chmod +x .claude/hooks/*.sh`
2. Verify `jq` is installed
3. Inspect hook execution with `claude --debug` or the transcript (Ctrl-R)
