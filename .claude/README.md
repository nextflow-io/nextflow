# Claude Code Hooks for Nextflow Development

This directory contains Claude Code hooks configured to improve the Nextflow development experience.

## Features

### 1. EditorConfig Enforcement (PostToolUse)
- **Trigger**: Immediately after editing any source file (`.groovy`, `.java`, `.gradle`, `.md`, `.txt`, `.yml`, `.yaml`, `.json`)
- **Action**: Applies editorconfig formatting rules using `eclint`
- **Files**: `hooks/format-editorconfig.py`

### 2. IntelliJ IDEA Formatter (PostToolUse)
- **Trigger**: Immediately after editing Groovy files
- **Action**: Applies IntelliJ IDEA code formatting to match project style
- **Files**: `hooks/format-idea.py`
- **Requirements**: IntelliJ IDEA installed (Community or Ultimate Edition)

### 3. Build Check (Stop + SubagentStop)
- **Trigger**: When Claude finishes responding or when a subagent completes
- **Action**: Runs `make compile` to verify code compiles without errors
- **Files**: `hooks/check-build.py`
- **Purpose**: Catch syntax and compilation errors immediately

### 4. Automatic Test Running (Stop)
- **Trigger**: When Claude finishes responding (main agent only)
- **Action**:
  - For source files: Runs corresponding test class (e.g., `CacheDB.groovy` → runs `CacheDBTest`)
  - For test files: Runs the specific test class
- **Files**: `hooks/run-tests.py`

## Hook Configuration

The hooks are configured in `.claude/settings.json` with:

### PostToolUse (runs after each edit)
- **format-editorconfig.py**: 30-second timeout
- **format-idea.py**: 60-second timeout

### Stop (runs when main agent finishes)
- **check-build.py**: 120-second timeout (runs first)
- **run-tests.py**: 300-second timeout (runs after build succeeds)

### SubagentStop (runs when subagent finishes)
- **check-build.py**: 120-second timeout

## Supported File Structure

The hooks understand Nextflow's module structure:

```
modules/
├── nextflow/src/main/groovy/nextflow/cache/CacheDB.groovy
├── nextflow/src/test/groovy/nextflow/cache/CacheDBTest.groovy
├── nf-commons/src/main/groovy/...
├── nf-lang/src/main/java/...
└── ...

plugins/
├── nf-amazon/src/main/nextflow/cloud/aws/...
├── nf-azure/src/main/nextflow/cloud/azure/...
└── ...
```

## Test Commands Generated

The hooks generate appropriate Gradle test commands:

- **Source file**: `modules/nextflow/src/main/groovy/nextflow/cache/CacheDB.groovy`
  - Runs: `./gradlew :nextflow:test --tests "*CacheDBTest"`

- **Test file**: `modules/nextflow/src/test/groovy/nextflow/cache/CacheDBTest.groovy`
  - Runs: `./gradlew :nextflow:test --tests "*CacheDBTest"`

- **Plugin file**: `plugins/nf-amazon/src/main/nextflow/cloud/aws/AwsPlugin.groovy`
  - Runs: `./gradlew :plugins:nf-amazon:test --tests "*AwsPluginTest"`

## Error Handling

- **EditorConfig failures**: Show warnings but don't block Claude
- **Test failures**: Provide detailed feedback to Claude for potential fixes
- **Missing tests**: Silently skip if no corresponding test exists
- **Timeouts**: Cancel long-running operations gracefully

## Dependencies

The hooks may automatically install:
- `eclint` via npm for editorconfig enforcement

Optional dependencies:
- IntelliJ IDEA (Community or Ultimate Edition) for Groovy formatting
  - Set `IDEA_SH` environment variable if not in standard location
  - Falls back gracefully if not found

## Hook Execution Flow

1. **During editing** (PostToolUse):
   ```
   ✓ EditorConfig formatting applied to CacheDB.groovy
   ✓ IDEA formatter: applied to CacheDB.groovy (5.2s)
   ```

2. **When finishing** (Stop):
   ```
   ✓ Build check passed (12.3s)
   ✓ Tests passed for CacheDB.groovy
   ```

3. **If errors occur**:
   ```
   Build check failed:
   error: cannot find symbol
     symbol:   variable foo
     location: class CacheDB
   ```

## Customization

You can modify the hooks by:
1. Editing the Python scripts in `hooks/`
2. Adjusting timeouts in `settings.json`
3. Adding or removing file extensions in the filter logic
4. Disabling specific hooks by removing them from `settings.json`

## Troubleshooting

If hooks aren't working:
1. Check that scripts are executable: `chmod +x .claude/hooks/*.py`
2. Verify Python 3 is available
3. Check Claude Code's debug output with `claude --debug`
4. Review hook execution in the transcript (Ctrl-R)
5. For IDEA formatter: verify IDEA installation or set `IDEA_SH` environment variable