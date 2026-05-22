# Native External Scripts for Nextflow Processes

- Authors: Edmund Miller, Nextflow maintainers
- Status: draft
- Date: 2026-05-22
- Tags: scripts, templates, modules, resource-bundles, typing, developer-experience

## Summary

Introduce a native-language `script:` execution mode for Nextflow processes where external Bash, Python, R, Julia, and related scripts are staged unchanged and receive task context through a typed runtime namespace, instead of being preprocessed as string templates. The feature should be designed as much for coding agents as for humans: agents should be able to inspect, edit, lint, and execute the same files that will run inside a Nextflow task.

## Problem Statement

Nextflow supports inline process scripts and external `template:` files, but `template:` files are not valid source files until Nextflow expands Groovy/GString-style placeholders. A Python, Bash, or R template can contain `$var` literals that are meaningful only during Nextflow template expansion, creating a file that is neither valid standalone source nor reliably understood by language-native tooling.

This creates recurring pain for both humans and coding agents:

- language servers, linters, type checkers, debuggers, and formatters cannot reason about the file as written;
- Bash templates require fragile escaping rules for values that should remain shell variables;
- large scripts are pushed back into triple-quoted process bodies because externalizing them loses tooling confidence;
- module authors have no clean, standard place for non-trivial scripts that are staged with the module and invoked with typed process inputs.

Agents can sometimes generate `template:` files quickly because templates are just text, but that is not the same as a maintainable source model. Template expansion hides the true runtime program from review, breaks language-native feedback loops, and makes it harder for humans to audit or reason about agent-authored changes. Native external scripts should make the agent-friendly path also be the human-readable path.

Snakemake's `script:` directive demonstrates a better shape: external scripts are real native-language files, and the workflow runtime injects a `snakemake` object or language-specific variables at execution time. Nextflow should borrow this pattern and improve it with Nextflow's typed process model.

## Goals or Decision Drivers

- Keep external script source files valid in their native language before, during, and after execution.
- Preserve IDE, linter, type checker, formatter, shellcheck, debugger, and unit-test workflows.
- Avoid template-time string rewriting for normal script context injection.
- Provide a consistent runtime namespace across Bash, Python, R, Julia, Rust, and future supported languages.
- Generate useful type stubs from the process signature so agents and IDEs can see typed inputs, outputs, params, and resources.
- Stage scripts through the same module/workflow resource mechanism as module `bin/` assets.
- Compose with `nextflow module run` so agents and developers can run a module directly, get fast feedback, and iterate on external script assets without writing a temporary wrapper workflow.
- Keep existing inline scripts and `template:` behavior compatible.

## Non-goals

- Remove or immediately deprecate existing inline process scripts.
- Remove the existing `template:` directive in the first implementation.
- Define a full package manager for script dependencies.
- Replace container, Conda, Wave, or module dependency semantics.
- Guarantee static typing for all dynamic channel shapes in the first release.

## Considered Options

- Keep `template:` as the only external script mechanism.
- Improve escaping and documentation for `template:`.
- Add native external scripts with runtime context injection.

## Pros and Cons of the Options

### Keep `template:` as the only external script mechanism

- Good, because it requires no new runtime model.
- Good, because existing pipelines continue unchanged.
- Bad, because template files remain invalid or misleading native-language source files.
- Bad, because it does not solve IDE, linting, testing, shell escaping, or agent-readability issues.

### Improve escaping and documentation for `template:`

- Good, because it reduces some user confusion.
- Good, because it is low implementation effort.
- Bad, because it treats the symptom rather than the model mismatch.
- Bad, because every language still has to coexist with Groovy string interpolation rules.

### Add native external scripts with runtime context injection

- Good, because script files are exactly the files that run.
- Good, because native-language tooling works without a Nextflow preprocessing step.
- Good, because script context can be typed and generated per process.
- Good, because it gives module authors and agents a clear default for non-trivial logic.
- Bad, because it requires runtime injection libraries/stubs and task sidecar metadata.
- Bad, because language-specific details must be designed carefully.

## Solution or decision outcome

Add a native external script mode for processes. Nextflow stages the script file unchanged into the task work directory and writes task context sidecars that language-specific shims expose as a `nextflow` namespace.

Example Bash script:

```bash
#!/usr/bin/env bash
set -euo pipefail

samtools sort -@ "${nextflow[cpus]}" \
  -o "${nextflow_output[bam]}" \
  "${nextflow_input[reads]}"
```

Example Python script:

```python
from nextflow.script import nextflow

samplesheet = nextflow.input["samplesheet"]
nextflow.output[0].parent.mkdir(parents=True, exist_ok=True)

pd.read_csv(samplesheet).to_parquet(nextflow.output[0])
```

The process declaration points to a native script asset rather than embedding a heredoc or preprocessing a template:

```groovy
process SORT_BAM {
  input:
  path reads

  output:
  path 'sorted.bam', emit: bam

  script:
  file 'scripts/sort_bam.sh'
}
```

The exact DSL spelling is open for design. If `script:` cannot be overloaded without ambiguity, alternatives include `scriptFile:`, `exec: file('scripts/sort_bam.sh')`, or `nativeScript:`. The architectural decision is the runtime model: stage unchanged source and inject context by sidecar, not by rewriting the source file.

## Rationale & discussion

### Snakemake precedent

Snakemake external scripts are regular files referenced from a rule. Its documentation describes access to a runtime `snakemake` object with input, output, params, wildcards, log, threads, resources, and config. Python scripts access values such as `snakemake.input[0]`; R scripts use an S4-style object; Bash scripts receive associative arrays such as `snakemake_input`, `snakemake_output`, `snakemake_params`, and `snakemake[threads]`.

The key design point is that the source file remains native-language source. Snakemake's own tests include real `test.py`, `test.R`, `test.sh`, and Julia/Rmd scripts under `tests/test_script/scripts`, exercising the runtime-injected namespace rather than template expansion. Snakemake also supports Rust scripts through `rust-script` and can expose a generated Rust `snakemake` instance from JSON task context, which is a useful precedent for typed compiled-language support.

Nextflow should adopt that separation of concerns:

1. the workflow engine resolves process inputs, outputs, params, resources, and task metadata;
2. the engine serializes that context into a task sidecar;
3. a small language-specific runtime exposes the context in native syntax;
4. the user script runs unchanged.

### Proposed task sidecars

For each task using native external scripts, Nextflow writes sidecars alongside `.command.sh`, for example:

- `.nextflow/task.json` — canonical task context: inputs, outputs, params, resources, environment, attempt, workDir, projectDir, moduleDir, container metadata, etc.;
- `.nextflow/bash-context.sh` — safely quoted Bash associative-array declarations;
- `.nextflow/nextflow_script.py` or an importable module path — Python loader and generated typing metadata;
- language-specific equivalents for R, Julia, and future languages.

The task wrapper sources or imports the appropriate shim before executing the staged script.

### Bash API

Bash should follow the Snakemake associative-array shape because it is simple and proven:

```bash
${nextflow_input[reads]}
${nextflow_input[0]}
${nextflow_output[bam]}
${nextflow_params[prefix]}
${nextflow[cpus]}
${nextflow[attempt]}
```

Guidelines:

- require Bash 4+ for associative arrays;
- quote all path-like values safely in the generated context file;
- represent nested values through JSON files or flattened keys rather than unsafe shell string interpolation;
- keep shell variable expansion in the user script native to Bash, not Nextflow.

### Python API

Python should provide an importable package and generated stubs:

```python
from nextflow.script import nextflow

nextflow.input["reads"]
nextflow.output.bam
nextflow.cpus
nextflow.params["prefix"]
```

A generated `.pyi` file can specialize the namespace for each process. If the process declares `path reads`, the stub can expose `nextflow.input["reads"]` as `pathlib.Path` or a list of `Path` when cardinality is known. `val` inputs can be typed from the DSL compiler where possible and fall back to `Any` where not.

### R, R Markdown, Julia, and Rust API

R and Julia should expose idiomatic objects while preserving the same conceptual fields:

- R: `nextflow@input[["reads"]]`, `nextflow@output[[1]]`, or a simpler list-like object if preferred;
- Julia: `nextflow.input["reads"]`, `nextflow.output[1]`, `nextflow.cpus`.

R Markdown (`.Rmd`) can be supported as an extension of the R runtime once the R context object is stable. The notebook/report use case is valuable, but it should not block the core external-script contract.

Rust should be listed as a future supported language. Snakemake executes Rust scripts via `rust-script` and requires `rust-script` plus OpenSSL and a C compiler toolchain in the rule environment; it then generates a typed Rust `snakemake` instance from JSON using `json_typegen`. Nextflow can follow the same broad model later: write the canonical task context sidecar, generate a typed Rust `Nextflow` struct from the process context, and run a valid `rust-script` source file unchanged when `rust-script` is available in the task environment.

The first implementation should prioritize Bash, Python, and R. Julia, Rust, R Markdown, and notebooks can be adopted later once the sidecar schema is stable.

### Future notebook integration

Jupyter notebooks should be an explicit future goal, not part of the initial native external script implementation. Snakemake's notebook integration is aimed at combining a modular workflow definition with small, focused notebooks for exploration and plotting rather than large monolithic notebooks. Nextflow can adopt the same product direction once the sidecar schema and Python runtime are mature:

- execute `.ipynb` assets as task scripts without rewriting notebook source;
- inject `nextflow` context into the notebook kernel through the same sidecar contract;
- preserve executed notebooks or rendered reports as declared outputs;
- keep notebooks small and module-local so they remain reviewable and reproducible.

### ResourceBundle integration

Native scripts should be staged through the same ResourceBundle mechanism proposed for module `bin/` and workflow assets:

- workflow scripts live in `scripts/` next to `bin/`;
- module scripts live inside the module directory, for example `modules/nf-core/bwa-align/scripts/sort_bam.py`;
- script-relative includes are resolvable from the script body through `scriptDir` or normal working-directory conventions;
- both local and registry-installed modules stage script assets consistently.

This makes external scripts a first-class module asset rather than an ad hoc path reference.

This is also the missing counterpart to module `bin/` support. `bin/` is useful for executable helper tools, but it does not give agents or maintainers a structured way to move process-specific logic out of heredocs while preserving typed process context. Native external scripts should cover that gap: helpers can remain in `bin/`, while task logic that depends on `input`, `output`, `params`, and resources can live in `scripts/` with an explicit runtime namespace.

### Agent feedback loop with `nextflow module run`

The module system ADR already proposes `nextflow module run scope/name` as a direct execution path for registry or local modules. Native external scripts should be designed to make that command a high-quality REPL-like feedback loop for agents and humans:

1. inspect a module's `main.nf`, `meta.yml`, `scripts/`, and `bin/` assets;
2. edit a native Bash/Python/R script using language-native tools;
3. run the module directly with representative inputs, for example `nextflow module run nf-core/bwa-align --reads reads.fq --reference genome.fa`;
4. inspect the task work directory, sidecar context, script outputs, and logs;
5. iterate without creating a throwaway wrapper workflow.

This should influence the implementation details:

- generated sidecars should be easy to find and read in the task work directory;
- errors should point back to the original script path and line number, not only to `.command.sh`;
- generated type stubs or context schemas should be available before task execution where possible;
- module-local `scripts/` and `bin/` assets should be staged consistently so an agent can reason about the complete runtime bundle.

### Relationship to `template:`

Native external scripts are a strict replacement for most uses of `template:` but do not need to break compatibility.

- Existing `template:` continues to perform legacy template expansion.
- Documentation should recommend native external scripts for non-trivial Bash/Python/R/Julia code.
- A future migration tool can convert simple templates to native scripts by replacing template variables with `nextflow.*` namespace references.
- Internally, legacy `template:` can eventually become a compatibility shim layered on top of ResourceBundle staging plus preprocessing.

### Agent and Co-Scientist impact

This feature gives recommendation systems and coding agents a positive default:

- keep small one-liners inline;
- move non-trivial process logic into native external scripts;
- lint and test those scripts directly;
- rely on generated stubs for typed I/O and autocomplete;
- use `nextflow module run` as the fast execution loop for module-local script changes;
- avoid heredoc/template escaping hazards.

That recommendation is more actionable than merely warning against large inline process bodies.

## Links

- Snakemake external script documentation: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts
- Snakemake script tests: https://github.com/snakemake/snakemake/tree/main/tests/test_script/scripts
- Related ADR: [Module System for Nextflow](20251114-module-system.md)
- Related ADR: [Typed Processes](20251017-typed-processes.md)
