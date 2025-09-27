# Claude Code Hooks for Nextflow Development

This directory contains Claude Code hooks configured to improve the Nextflow development experience.

## Features

### 1. EditorConfig Enforcement
- **Trigger**: After editing any source file (`.groovy`, `.java`, `.gradle`, `.md`, `.txt`, `.yml`, `.yaml`, `.json`)
- **Action**: Applies editorconfig formatting rules using `eclint`
- **Files**: `hooks/format-editorconfig.py`

### 2. Automatic Test Running
- **Trigger**: After editing source files or test files in modules or plugins
- **Action**:
  - For source files: Runs corresponding test class (e.g., `CacheDB.groovy` → runs `CacheDBTest`)
  - For test files: Runs the specific test class
- **Files**: `hooks/run-tests.py`

## Hook Configuration

The hooks are configured in `.claude/settings.json` with:
- **30-second timeout** for editorconfig formatting
- **5-minute timeout** for test execution
- **Smart file filtering** to only process relevant files
- **Parallel execution** of both hooks after file edits

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

## Customization

You can modify the hooks by:
1. Editing the Python scripts in `hooks/`
2. Adjusting timeouts in `settings.json`
3. Adding or removing file extensions in the filter logic

## Troubleshooting

If hooks aren't working:
1. Check that scripts are executable: `chmod +x .claude/hooks/*.py`
2. Verify Python 3 is available
3. Check Claude Code's debug output with `claude --debug`
4. Review hook execution in the transcript (Ctrl-R)

## Example Output

When editing a file, you'll see:
```
✓ EditorConfig formatting applied to CacheDB.groovy
✓ Tests passed for CacheDB.groovy
BUILD SUCCESSFUL in 2s
```

If tests fail:
```
Tests failed for CacheDB.groovy:
Error output:
CacheDBTest > testCacheCreation FAILED
    AssertionError: Expected true but was false
```