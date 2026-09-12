# Process error section

- Authors: Ben Sherman
- Status: proposed
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-06-15
- Tags: lang, static-types, processes, error-handling

## Summary

Add an `error:` section to typed processes, mirroring the `output:` section, that lets a process emit a *domain error* as a value instead of aborting the run.

## Problem Statement

When a task fails, Nextflow handles it using the `errorStrategy` directive, which provides a few standard behaviors such as retrying the task, ignoring the error, or terminating the pipeline. In the case of termination, Nextflow reports the error with a standard format that accounts for all of the various ways in which a task can fail.

However, a task can fail for two fundamentally different reasons:

- **Execution errors** -- the infrastructure failed: out of memory, node lost, spot reclaim, submit failure. These are transient or environmental and are typically handled by retries.

- **Domain errors** -- the tool ran correctly but legitimately failed on the *data*: an input sample was too low quality to align, no variants were called, a file was malformed. These are an expected outcome for some inputs, not a fault of the pipeline or the infrastructure.

There are currently two ways to handle domain errors manually:

- **Use the `ignore` error strategy.** This allows the pipeline to complete, but it does not provide an easy way to manage failed inputs. The developer must join the input channel to the output channel and filter out values for which the output is missing, and there is no way for a failed task to return anything (e.g. an error log).

- **Catch domain errors in the process script.** The developer can add logic in the process script to catch certain error conditions, so that they can still output something. This requires clobbering the process outputs to handle both success and error conditions (e.g. making all outputs optional, providing "fake" outputs on failure) and filtering output channels to separate successful and failed tasks.

The first approach is a quick fix, not a real error handling solution. The second approach is nearly a complete solution, but lacks a key ingredient -- the ability to emit a separate "error" output for domain errors.

## Goals

- Allow a process to emit a domain error as a value, so the workflow can handle it with normal dataflow logic.

- Maintain backwards compatibility with existing code -- domain errors should be a per-process, opt-in feature.

- Continue to use `errorStrategy` for execution errors.

## Non-goals

- Support for legacy processes. The `error:` section is introduced for typed processes only; legacy support may follow later.

- New syntax for *classifying* errors (no exit-code lists, no guard expressions).

## Solution

Introduce an **`error:` section** for typed processes, with the same syntax as the `output:` section. A domain error is detected structurally -- **a task that succeeds (exit 0) but does not fulfill its declared `output:`** -- and is emitted as an error value instead of triggering the error strategy.

## Core Capabilities

### Domain errors vs execution errors

When a task completes, the following strategy is used to detect domain errors vs execution errors:

1. **Exit code != 0** → execution error → `errorStrategy`.
2. **Exit 0, `output:` fulfilled** → emit `output:` (normal success). The output path wins even if error artifacts also happen to be present.
3. **Exit 0, `output:` not fulfilled, `error:` fulfilled** → domain error → emit `error:`, the workflow continues, `errorStrategy` is *not* triggered.
4. **Exit 0, `output:` not fulfilled, `error:` not fulfilled** → execution error → `errorStrategy`.

The output section is considered **not fulfilled** if any required output files (`file()`, `files()`) or environment variables (`env()`) are missing. Other errors such as a missing variable, missing `stdout()`, or missing `eval()` are not treated as domain errors because they usually indicate a malformed pipeline.

Therefore, the only way to trigger a domain error is to exit 0, ensure that the normal outputs are missing, and ensure that the error outputs are present. The pipeline developer is responsible for writing the process in this way, for example:

```groovy
nextflow.enable.types = true

process ALIGN {
    input:
    record(id: String, reads: Path)
    index: Path

    output:
    record(id: id, bam: file('aligned.bam'))

    error:
    record(id: id, log: file('aligner.log'))

    script:
    """
    aligner ${reads} ${index} > aligned.bam 2> aligner.log || {
        rm -f aligned.bam
        exit 0
    }
    """
}
```

The `error:` section is optional. However, the `output:` section is required when `error:` is defined.

### Process call semantics

A process that declares both `output:` and `error:` has return type `Tuple<V, E>`, where `V` is the output type and `E` is the error type. The caller should destructure the tuple to access output and error separately:

```groovy
ch_smaples = channel.of( /* ... */ )
(ch_aligned, ch_failed) = ALIGN(ch_samples)

ch_aligned.view { r -> "Aligned sample ${r.id}: ${r.bam}" }
ch_failed.view { r -> "Failed to align sample ${r.id}: ${r.log}" }
```

When a process with both `output:` and `error:` returns a dataflow channel, it emits each task result to either the output channel or error channel depending on whether the result is a domain error. Thus, the total number of output values and error values is always equal to the number of inputs.

When a process returns a dataflow value, the output and error values are each bound to either the task result or `null`, depending on whether the result is a domain error. Thus, the dataflow values for output and error are always bound to a value (no "empty" dataflow value).

### Execution, caching, and lineage

A task that emits an `error:` is treated as a **successful, cacheable task** (exit 0). On a resumed run, the domain error is cached and the task is not re-executed. The cache entry does not need to be modified, since the domain error is re-derived from the task directory (outputs missing, error outputs present).

The `TaskOutput` lineage record should be extended with an `error` field that mirrors the existing `output` field. These fields should be mutually exclusive.

When a domain error occurs, `topic:` emissions and `publishDir` are skipped.

## Alternatives

### Triggering domain errors with exit codes

One alternative is to trigger domain errors by returning certain exit codes. Many command-line tools use exit codes for this very purpose, e.g. to distinguish an invalid input from an out-of-memory error. In fact, earlier versions of Nextflow had a `validExitStatus` process directive for this very purpose.

However, this approach does not work in general:

- Different tools use different exit code conventions.
- There is no way to know which command in a script returned the exit code.

The `validExitStatus` directive was ultimately removed for these same reasons. While it seems intuitive to simply rely on exit codes, this interface is not rich enough to classify domain vs. execution errors.

Instead, pipeline developers must write process scripts in a way that triggers domain errors when desired. This approach is more verbose, but it seems to be the only one that works across all possible tools and environments.

### Triggering domain errors via `emit` error strategy

Sometimes, it is useful to treat execution errors as domain errors for practical reasons. For example, given a task that runs out of memory even after multiple attempts with additional memory, the user might want to treat this task as a "lost cause" so that the rest of the pipeline can proceed.

This can be achieved by adding an `emit` error strategy which simply emits the task failure as a domain error using the `error:` section:

```groovy
nextflow.enable.types = true

process ALIGN {
    memory { 8.GB * task.attempt }
    errorStrategy { task.attempt < 3 ? 'retry' : 'emit' }

    input:
    record(id: String, reads: Path)
    index: Path

    output:
    record(id: id, bam: file('aligned.bam'))

    error:
    record(id: id, log: file('aligner.log'))

    script:
    """
    aligner ${reads} ${index} > aligned.bam 2> aligner.log || {
        rm -f aligned.bam
        exit 0
    }
    """
}
```

However, the `error:` section is not reliable in the event of an execution error, since the task could have failed before the error outputs were created. In that case, the error strategy would have to fallback to a default strategy, likely `terminate` or `finish`.

As a result, it is unclear whether an `emit` strategy would actually be useful. It remains a possibility for future investigation.

### Handling domain errors with try/catch/throw

Many languages, including Java and Python, provide the ability to *throw* or *raise* an error up the call stack. Any upstream caller can *catch* this error and handle it; otherwise it is handled by the runtime. This approach is flexible because it allows errors to be propagated through multiple levels of indirection with minimal ceremony.

Nextflow inherits the try-catch-throw syntax from Java and Groovy, mainly for compatibility with existing Nextflow code and Java/Groovy libraries that can throw exceptions. However, Nextflow's dataflow programming model fits much better with errors-as-values because it provides a clear flow of data that works across any level of concurrency.

## Links

- Community issues: [#725](https://github.com/nextflow-io/nextflow/issues/725), [#903](https://github.com/nextflow-io/nextflow/issues/903)
- Related: [Typed processes](20251017-typed-processes.md)
- Related: [Typed workflows](20260310-typed-workflows.md)
