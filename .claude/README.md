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

### 2. Build check (Stop + SubagentStop)
- **Script**: `hooks/check-build.sh`
- **Trigger**: when Claude (or a subagent) finishes responding
- **Action**: runs `make compile` to catch syntax/compilation errors early
- **Failure mode**: blocking — feeds a parsed error summary back to Claude so it
  can fix the breakage before yielding

### 3. Test runner (PostToolUse)
- **Script**: `hooks/run-tests.sh`
- **Trigger**: after an `Edit`/`Write`/`MultiEdit` on a `.groovy`/`.java` file under
  `modules/` or `plugins/`
- **Action**: maps the edited file to its test class and runs it, e.g.
  `modules/nextflow/.../cache/CacheDB.groovy` → `./gradlew :nextflow:test --tests "*CacheDBTest"`
  (test files run themselves)
- **Failure mode**: blocking — feeds the test-failure summary back to Claude
- **Scoping**: the script itself decides what is testable (`.groovy`/`.java` under
  `modules/` or `plugins/`, mapped to a gradle module). This covers all three edit tools
  without gaps, so no declarative `if` is used here. Runs per edit, so the timeout is
  generous (300s); narrow the matcher if per-edit test runs are too heavy.

### 4. IntelliJ IDEA formatter (PostToolUse, async)
- **Script**: `hooks/format-idea.sh`
- **Trigger**: after an `Edit`/`Write`/`MultiEdit` on a `.groovy` file
- **Action**: runs IntelliJ IDEA's headless formatter to match project code style
- **Dependencies**: IntelliJ IDEA (Community or Ultimate). Located via the `IDEA_SH`
  environment variable or common macOS/Linux/JetBrains-Toolbox install paths
- **`async: true`**: runs in the background; non-blocking — skips silently if IDEA is not
  installed and only warns (never blocks) if the formatter errors

### 5. Integration test runner (PostToolUse, asyncRewake)
- **Script**: `hooks/run-integration-test.sh`
- **Trigger**: editing an integration case `tests/<name>.nf` (scoped declaratively via
  `if: "Edit(tests/*.nf)"` / `Write(...)` / `MultiEdit(...)` — one entry per edit tool,
  since `if` patterns are tool-qualified)
- **Action**: runs just that one case through the existing runner
  (`tests/checks/run.sh ../<name>.nf`) rather than the full ~115-case suite. Uses the
  dev launcher (`launch.sh`) and enables Docker only when it is available.
- **`asyncRewake: true`**: runs in the background (never blocks the conversation) and, on
  failure (exit code 2), wakes Claude with the failure output so it can react.
- **Requirements**: a built dev launcher (`make compile`) and, for Docker-only cases,
  Docker. The runner writes scratch/output dirs under `tests/checks/`; `tests/cleanup.sh`
  clears them.

## Enabling / disabling

Hooks are configured in `.claude/settings.json`. Disable one by removing its entry;
adjust timeouts there as needed.

## Troubleshooting

1. Ensure scripts are executable: `chmod +x .claude/hooks/*.sh`
2. Verify `jq` is installed
3. Inspect hook execution with `claude --debug` or the transcript (Ctrl-R)
