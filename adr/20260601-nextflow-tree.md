# Agent-friendly project outline with `nextflow tree`

- Authors: Edmund Miller, OpenAI
- Status: draft
- Date: 2026-06-01
- Tags: cli, lang, parser, tooling, agents

## Summary

Add a `nextflow tree` command that prints a concise structural outline of a Nextflow project for coding agents and humans. The command shows workflows, subworkflows, modules/processes, and their call relationships, while omitting executable syntax and low-level implementation details.

## Problem Statement

Coding agents often need to understand the shape of a Nextflow project before editing it. Today they must read many `.nf` files, includes, process bodies, channel expressions, and config files to infer a simple question: “what workflows call what?” This is slow, noisy, and error-prone.

Nextflow already has the semantic information needed to answer this question: the strict parser exposes workflow, process, and include nodes, and include resolution can map aliases to target components. A dedicated CLI command should present that information in a stable, compact, agent-friendly format.

## Goals

- Show the high-level project call graph rooted at one or more workflows.
- Include workflows, subworkflows, included components, process/module calls, aliases, and unresolved references.
- Omit channel expressions, process scripts, directives, params, `take:`, `main:`, `emit:`, and other syntax that does not help understand structure.
- Be useful by default on incomplete or partially broken projects.
- Provide a stable text format for agents and optional machine-readable output.
- Reuse Nextflow parser and include-resolution semantics rather than inventing a separate language model.

## Non-goals

- Execute the pipeline.
- Validate dataflow correctness or types beyond what parsing and include resolution already provide.
- Replace `nextflow inspect`, linting, language-server outline views, or documentation generators.
- Show all dependencies, operators, channels, params, config profiles, containers, or resources.
- Require projects to use a specific module layout beyond existing Nextflow include rules.

## Decision

Implement `nextflow tree [SCRIPT]` as a best-effort structural outline command.

Default invocation:

```bash
nextflow tree
```

is equivalent to inspecting `main.nf` in the current project and printing the entry workflow plus any statically resolvable subworkflow and process/module calls.

Example output:

```text
workflow main
  subworkflow preprocess
    module FASTQC
    module TRIM_GALORE
  subworkflow align
    module BWA_MEM
    module SAMTOOLS_SORT
  module MULTIQC
```

The command should use the strict parser / `nf-lang` model as the source of truth. In particular, it should build on script AST concepts such as `WorkflowNode`, `ProcessNode`, and `IncludeNode`, and include-resolution behavior such as `ResolveIncludeVisitor`, which already knows how to resolve local includes, directory includes, remote modules, plugin includes, aliases, and missing included names.

## User interface

```text
nextflow tree [SCRIPT]

Options:
  --workflow <name>       Show tree rooted at a named workflow
  --all                   Show all workflows, not only the entry workflow
  --depth <n>             Limit recursive expansion depth
  --files                 Include source file paths
  --aliases               Show alias/original-name details
  --unresolved            Include unresolved calls/includes in the tree
  --no-expand-includes    Show included components but do not recursively inspect their files
  --json                  Emit machine-readable JSON
  --dot                   Emit Graphviz DOT
  --strict                Exit non-zero on parse or resolution errors
```

Default behavior should be forgiving: print the best recovered tree and append warnings. `--strict` changes warnings that affect correctness into errors.

## Text output model

The default text output has three node labels:

```text
workflow <name>
subworkflow <name>
module <name>
```

Rules:

- The unnamed entry workflow is printed as `workflow main`.
- Named workflows at the root are printed as `workflow <name>`.
- Workflow calls nested inside another workflow are printed as `subworkflow <name>`.
- Process calls and included process components are printed as `module <name>`.
- Alias names are the primary visible names.
- Indentation represents call nesting.
- Repeated calls are preserved when they appear in distinct structural positions.
- Warnings are printed after the tree, not inline, unless `--unresolved` is enabled.

## Details intentionally omitted

The command should not display:

- full Nextflow syntax;
- `take:`, `main:`, `emit:`, `publish:`, `workflow.onComplete`, or `workflow.onError` sections;
- channel construction and operators;
- process inputs, outputs, directives, and script bodies;
- params and config values;
- assignments and intermediate variables, except when needed internally for name resolution;
- function bodies;
- comments;
- container, resource, executor, or profile details.

For example:

```nextflow
workflow {
  reads = Channel.fromPath(params.reads)
  FASTQC(reads)
  TRIM_GALORE(reads)
  MULTIQC(FASTQC.out.zip)
}
```

prints as:

```text
workflow main
  module FASTQC
  module TRIM_GALORE
  module MULTIQC
```

## Resolution behavior

### Includes and aliases

Given:

```nextflow
include { FASTQC; MULTIQC as QC_REPORT } from './modules/qc'
```

Default output:

```text
workflow main
  module FASTQC
  module QC_REPORT
```

With `--aliases`:

```text
workflow main
  module FASTQC
  module QC_REPORT  # MULTIQC from ./modules/qc
```

### Subworkflows

If a workflow calls another workflow, expand it recursively:

```text
workflow main
  subworkflow rnaseq
    module STAR_ALIGN
    module FEATURECOUNTS
  module MULTIQC
```

If recursion would cycle, show the repeated subworkflow once and mark it:

```text
workflow main
  subworkflow loop
    subworkflow loop  # recursive
```

### Duplicate modules

If the same module is called more than once, preserve repeated calls:

```text
workflow main
  module ALIGN
  module ALIGN
```

If aliases distinguish repeated imports, show the aliases:

```text
workflow main
  module ALIGN_TUMOR
  module ALIGN_NORMAL
```

With `--aliases`:

```text
workflow main
  module ALIGN_TUMOR   # ALIGN from ./modules/align
  module ALIGN_NORMAL  # ALIGN from ./modules/align
```

If a name is ambiguous, keep the tree readable and warn:

```text
workflow main
  ambiguous module FASTQC

warnings:
  FASTQC is defined or imported more than once; use --aliases or --files for details
```

### Unresolved includes and calls

Default output should include warnings:

```text
workflow main
  module FASTQC

warnings:
  unresolved include: ./modules/reporting
  unresolved call: ALIGN
```

With `--unresolved`, unresolved nodes are included inline:

```text
workflow main
  module FASTQC
  unresolved module ALIGN
  unresolved include ./modules/reporting
```

Dynamic calls that cannot be resolved statically should not be hidden:

```text
workflow main
  unresolved call <dynamic>
```

### Partial parse failures

Default mode should recover whatever structure is available:

```text
workflow main
  module FASTQC
  unresolved subworkflow report

warnings:
  partial parse failure in ./subworkflows/report.nf; showing recovered outline only
```

With `--strict`, the command exits non-zero:

```text
error: failed to parse ./subworkflows/report.nf
```

## JSON output

`--json` should expose a stable tree with warnings:

```json
{
  "workflows": [
    {
      "name": "main",
      "file": "main.nf",
      "children": [
        {
          "type": "subworkflow",
          "name": "preprocess",
          "file": "./subworkflows/preprocess.nf",
          "children": [
            { "type": "module", "name": "FASTQC", "file": "./modules/fastqc/main.nf" },
            { "type": "module", "name": "TRIM_GALORE", "file": "./modules/trim_galore/main.nf" }
          ]
        },
        { "type": "module", "name": "MULTIQC", "file": "./modules/multiqc/main.nf" }
      ]
    }
  ],
  "warnings": []
}
```

JSON fields should be additive over time. Initial fields:

- `type`: `workflow`, `subworkflow`, `module`, `unresolved`, or `ambiguous`;
- `name`;
- `originalName`, when aliased and requested or useful;
- `file`, when known;
- `line`, when known;
- `children`;
- `warnings`.

## Implementation sketch

1. Parse the requested script and reachable included scripts with `nf-lang`.
2. Run include resolution using the same rules as the strict parser:
   - relative file includes;
   - directory includes resolving to `main.nf`;
   - omitted `.nf` suffixes;
   - remote module references;
   - plugin includes;
   - aliases.
3. Build a symbol table per script:
   - workflows;
   - processes;
   - included components;
   - alias to target mappings;
   - source locations.
4. Visit workflow bodies and collect only top-level callable invocations.
5. Classify each invocation:
   - target is workflow: `subworkflow`;
   - target is process: `module`;
   - target is missing: unresolved;
   - target has multiple candidates: ambiguous;
   - target is dynamic: unresolved dynamic call.
6. Recursively expand subworkflows until `--depth`, cycle detection, or unresolved boundary.
7. Render text, JSON, or DOT.

This can start as a static parser feature. Later, it can share more code with the language server outline/call-hierarchy implementation if those APIs become public.

## Example project shapes

### Single-file pipeline

```text
workflow main
  module FASTQC
  module MULTIQC
```

### Nested subworkflows

```text
workflow main
  subworkflow preprocess
    module FASTQC
    module TRIM_GALORE
  subworkflow quantify
    module SALMON_INDEX
    module SALMON_QUANT
  module MULTIQC
```

### Multiple workflows

```text
workflow main
  subworkflow rnaseq
    module STAR_ALIGN
    module FEATURECOUNTS

workflow test
  module MAKE_TEST_READS
  subworkflow rnaseq
    module STAR_ALIGN
    module FEATURECOUNTS
```

### Include boundary without expansion

With `--no-expand-includes`:

```text
workflow main
  subworkflow preprocess  # from ./subworkflows/preprocess
  module MULTIQC          # from ./modules/multiqc
```

### Broken but useful project

```text
workflow main
  module FASTQC
  unresolved module ALIGN
  unresolved subworkflow report

warnings:
  missing include target: ./modules/align
  failed to parse ./subworkflows/report.nf
```

## Open questions

- Should `module` mean only included process components, or all process calls including local processes? This ADR recommends using `module` for all process calls because it is more useful to agents.
- Should the default root include all workflows or only the executable entry workflow? This ADR recommends entry workflow only, with `--all` for complete listings.
- Should `nextflow tree` live in core CLI or a plugin first? Core CLI is preferred because it depends on parser semantics and should be available to agents without project-specific setup.
