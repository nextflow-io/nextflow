# Runtime config origins

- Authors: Phil Ewels
- Status: draft
- Date: 2026-11-02
- Tags: config, cli, diagnostics

## Summary

Record where every resolved config value came from while the config is built: file and line, profile, selector, command-line flag, params file, environment variable, or default. Show it where users need it: the failed task report, `nextflow config -show-origin`, config warnings, and lineage.

## Problem statement

A setting in a Nextflow run can be set in many places, and the runtime applies them in a fixed order:

| Source | Example | Applied in |
| --- | --- | --- |
| Home config | `$NXF_HOME/config` | `ConfigCmdAdapter.resolveConfigFiles` |
| Pipeline config | `<projectDir>/nextflow.config` | `ConfigCmdAdapter.resolveConfigFiles` |
| Launch dir config | `./nextflow.config` | `ConfigCmdAdapter.resolveConfigFiles` |
| User config | `-c custom.config` | `ConfigCmdAdapter.resolveConfigFiles` |
| Includes, local or remote | `includeConfig 'https://raw.githubusercontent.com/nf-core/configs/...'` | `ConfigDsl.includeConfig` |
| Profiles, in `-profile` order | `profiles { uppmax { ... } }` | `ProfilesDsl.apply` |
| Command-line options | `-process.memory`, `-with-docker`, `-work-dir` | `ConfigCmdAdapter.configRunOptions` |
| Params | `params` block default, config, `-params-file`, `--name` | `CmdRun.parsedParams`, `ConfigDsl.withCliOverride`, `ParamsDsl` |
| Environment | `NXF_WORK` as the `workDir` default | `ConfigCmdAdapter.configRunOptions` |
| Process layers | process body, `withLabel`, `withName` (by base, simple, then fully qualified name), process scope defaults | `ProcessConfigBuilder.applyConfig` |
| Per task | dynamic values such as `{ 6.GB * task.attempt }`, capped by `resourceLimits` | `TaskConfig.getMemory` and friends |

The runtime merges all of this into one `ConfigObject` and keeps no record of where each value came from. That has three costs.

**Failed tasks don't say what to change.** "How do I give process X more memory?" is one of the most common support questions in nf-core. When a task fails, `TaskErrorFormatter.formatTaskError` prints the command, exit status, stdout, stderr, work dir and container. It doesn't print the resources the task actually had, or which setting produced them. To answer that question, the user has to reproduce the precedence rules by hand across `conf/base.config`, `conf/modules.config`, an institutional profile fetched from a URL, and their own `-c` file.

**`nextflow config` shows the result, not the reason.** It prints the merged config but not which file, line or flag set each value.

**Config warnings don't say where.** `Unrecognized config option 'process.foo'` (`ConfigValidator.validateOption`) doesn't name the file or line.

The VS Code extension 1.8.0 added a config preview for a process (language server `ConfigPreviewProvider`). For each directive it lists the layers that set it, with file, line and profile, and marks which one wins. That answers the question in the editor, but only statically: it renders source text without evaluating it, and it only sees config files in the workspace. It can't see `-c` files, remote includes, command-line options, params files, environment variables, values computed per task, or `resourceLimits` capping. Those are the cases users hit at runtime. It also re-implements selector matching and precedence (its comments point at `REPEATABLE_DIRECTIVES` in `ProcessConfigBuilder`), so the two implementations can drift.

## Goals

- Every resolved config option can report the ordered list of assignments that set it, with the one in effect marked.
- Every task directive can report the layer that set it (process scope, process body, `withLabel:X`, `withName:Y`, command line). Resource directives also report whether `resourceLimits` capped them, and where that limit was set.
- Show this in the failed task report, `nextflow config` and config warnings, and record it in lineage.
- One implementation of selector matching, shared by the runtime and the language server.
- No behavior change. Values, precedence and task hashes are unaffected, and recording origins can't fail a run.
- Negligible cost. No extra file reads, and work done per process rather than per task.

## Non-goals

- The v1 config parser (`NXF_SYNTAX_PARSER=v1`). Origins are reported as unknown.
- Tracing through expressions. `memory = params.mem` reports the line of that assignment, not where `params.mem` was set, although the param's own origin is available separately.
- Changing precedence rules.
- Platform and language server UI changes. This ADR defines the data and a JSON output for them to consume.
- Directives in lineage task records. `TaskRun` doesn't record directives today, and adding them is a separate change.

## Considered options

### A. Side table of origins, recorded at assignment

Each compiled config assignment carries its source line. `ConfigDsl` appends an origin to a side table keyed by the option path. Config values themselves are unchanged.

- Good, because there is one place that sees every config value. The config compiler turns every assignment into `assign(names, value)` (`ConfigToGroovyVisitor.transformConfigAssign`, `ConfigToGroovyVisitor.java:91`), and every such call goes through `ConfigDsl.assign` (`ConfigDsl.groovy:150`).
- Good, because nothing that reads config values has to change, including plugins.
- Good, because it captures evaluated values, remote includes, profiles and command-line options.
- Bad, because the side table must be merged everywhere the config is merged: config files, includes, profiles and direct writes from command-line options.
- Bad, because it only works with the v2 parser.

### B. Wrap each value in an origin-carrying type

Similar to Spring Boot's `OriginTrackedValue`.

- Good, because a value and its origin can't get out of sync.
- Bad, because every reader of `session.config` and every plugin must unwrap values, which breaks the plugin API.
- Bad, because Groovy truthiness and equality on wrapped values would become a source of subtle bugs.

### C. Re-derive origins statically on demand

Re-parse the config files with nf-lang when asked, as the language server does.

- Good, because the runtime doesn't change.
- Bad, because it can't evaluate values, and misses command-line options, params files and environment variables. That is exactly the gap this ADR is about.
- Bad, because it adds a second copy of the precedence rules to the runtime.

### D. Log every assignment at trace level

- Good, because it's trivial.
- Bad, because it can't be queried or used in error messages, and it floods the log.

## Solution

Adopt option A. Origins are recorded while the config is built, kept on the `Session` next to the config, resolved per process by `ProcessConfigBuilder`, and shown in four places: the failed task report, `nextflow config -show-origin`, config warnings and lineage.

## Rationale and discussion

### Data model

```groovy
class ConfigOrigin {
    enum Kind { FILE, CLI, PARAMS_FILE, ENV, SCRIPT, DEFAULT, UNKNOWN }
    Kind kind
    String source     // file URI, flag (`--outdir`, `-process.memory`), env var name, or script path
    int line          // 1-based, or -1 when not applicable
    String profile    // the profile block the assignment was in, or null
}

class ConfigOrigins {
    List<ConfigOrigin> history(List<String> path)   // every assignment, in the order applied
    ConfigOrigin effective(List<String> path)       // the last one; UNKNOWN rather than null
    void record(List<String> path, ConfigOrigin origin)
    void merge(ConfigOrigins other, List<String> prefix)
}
```

Paths are lists rather than dotted strings, because selector keys can contain dots and colons (`withName:'.*:STAR_ALIGN'`).

The table keeps the full history, not only the value in effect. That lets every consumer show what was overridden, as the VS Code preview does.

A map assignment such as `resourceLimits = [memory: 72.GB]` gets one origin, and its keys inherit it.

### Recording

| Where values enter | Change |
| --- | --- |
| Config assignments | `ConfigToGroovyVisitor.transformConfigAssign` passes the node's line: `assign(names, value, line)`. `ConfigDsl.assign` and the block DSLs record `FILE(configPath, line, currentProfile)`. |
| Profiles | `ProfilesDsl.apply` sets the current profile name on the DSL while it replays each profile block, and clears it afterwards. |
| `includeConfig` | The nested parser's origins are merged under the include's scope prefix, next to the existing `deepMerge`. The source is `includePath.toUriString()`, so a remote include shows its URL. |
| Config files | `ConfigBuilder.merge0` merges each file's origins after `result.merge(config)`. |
| `-cluster.*`, `-executor.*` | These are parsed as config text. They are given a source label naming the flag instead of a path. |
| Other command-line options | `ConfigCmdAdapter.configRunOptions` records `CLI(<flag>)` for each key it writes, and `ENV(NXF_WORK)` and similar for environment defaults. |
| Params | `CmdRun.parsedParams` tags each param as `PARAMS_FILE` or `CLI`. `ConfigDsl.withCliOverride` records that origin when it overrides a config assignment. The script `params` block records `SCRIPT` defaults. |

Config files are compiled on every run, so changing the generated `assign` call only affects the DSL classes in this repository. A two-argument overload stays, with origin `UNKNOWN`.

### Per process and per task

`ProcessConfigBuilder.applyConfig` already applies the layers in order: process body, matching `withLabel` selectors, matching `withName` selectors, then process scope defaults for anything still unset. `applyConfigSettings` gains the path of the settings it applies (for example `['process', 'withLabel:process_high']`), so for each directive it writes it can record a `DirectiveOrigin(layer, rule, ConfigOrigin)` on the `ProcessConfig`. Directives in the process body are recorded as `SCRIPT` with the process name.

`TaskConfig` gains `getOrigin(String directive)`. For `cpus`, `memory`, `time` and `disk` it also reports when `resourceLimits` capped the value, and the origin of that limit.

A dynamic directive keeps the origin of its closure. The value shown is the one evaluated for the current attempt.

### One selector implementation

`ConfigSelectorResolver` in `modules/nextflow/src/main/groovy/nextflow/script/dsl` depends only on `java.util.regex`. Move it, together with the list of repeatable directives, to nf-lang, so the language server's config preview uses the same matching and ordering as the runtime.

### Where origins are shown

**1. Failed task report.** `TaskErrorFormatter.formatTaskError` adds a `Task resources:` block before the tip, listing only the directives that are set:

```
Task resources:
  cpus   = 12     withLabel:process_high  conf/base.config:41
  memory = 72 GB  withLabel:process_high  conf/base.config:42  (attempt 2)
           capped by resourceLimits.memory  profile uppmax  https://raw.githubusercontent.com/nf-core/configs/master/conf/uppmax.config:18
  time   = 16h    withLabel:process_high  conf/base.config:43
```

The same text goes into `WorkflowMetadata.errorReport`, which Seqera Platform receives, and is passed to `TaskTipProvider` extensions, so a tip can name the selector to override.

**2. `nextflow config -show-origin`.** The canonical and flat outputs append the origin in effect as a comment:

```
process {
   resourceLimits = [cpus: 16, memory: 72.GB, time: 24.h]   // profile uppmax, https://raw.githubusercontent.com/nf-core/configs/master/conf/uppmax.config:18
   withLabel:process_high {
      cpus = { 12 * task.attempt }         // conf/base.config:41
      memory = { 36.GB * task.attempt }    // conf/base.config:42
   }
   withName:STAR_ALIGN {
      memory = 100.GB                      // custom.config:3 (-c)
   }
}
params {
   input = 'samplesheet.csv'               // params.yaml (-params-file)
   outdir = 'results'                      // --outdir
}
```

`-o json -show-origin` prints a list of `{path, value, origins}` entries. This is the machine-readable form for the language server and Platform. `-value <option> -show-origin` prints the full history for one option.

**3. Config warnings.** `ConfigValidator`'s `Unrecognized config option` warning includes the file and line.

**4. Lineage.** `WorkflowRun` already records `config`. It gains `configOrigins` in the JSON form above, so `nextflow lineage diff` between two runs shows where each side's value came from.

### Cost

One small object per assignment. Typical nf-core pipelines have hundreds to a few thousand config assignments. Directive origins are resolved once per process, not per task. No files are read beyond those already parsed.

### Risks

- **Missing origins.** A new code path that writes into the config without recording an origin shows `UNKNOWN`. To catch this, `effective()` never returns null, and a test asserts that every option in the integration test configs has a known origin.
- **Leaking credentials.** No config values are printed anywhere new except the resource block, which only holds resource values. Origins are paths and URLs, and URLs are printed without user info or query string, because a private include URL can carry a token.
- **Longer error reports.** The resource block is at most five lines.

## Open questions

- Should `Task resources:` appear on every failure, or only when the exit status or executor reason points at resources (137, 140, 143, out of memory, timeout)? The proposal is every failure: it's a few lines, and resource misconfiguration shows up under many exit codes.
- A per-process view on the command line (`nextflow config -process STAR_ALIGN`) needs the process labels, which are declared in the script. Should `nextflow config` compile the script for this, as `-preview` does, or should the per-process view stay in the language server?
- Should `SCRIPT` origins for process body directives include line numbers in this phase? That needs the script compiler to pass positions for directive calls, the same change as for config assignments.

## Links

- [VS Code extension: config preview for a process](https://github.com/nextflow-io/vscode-language-nextflow/issues/240)
- Depends on the v2 config parser from [Strict syntax parser](20250508-strict-syntax-parser.md)
- Params precedence from [Workflow params](20250825-workflow-params.md)
- Lineage records from [Data lineage](20250508-data-lineage.md)
- Prior art: `git config --show-origin`, Spring Boot configuration origin tracking
