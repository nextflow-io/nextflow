
# Package `nextflow.executor`

The `nextflow.executor` package defines the executor interface and implements several built-in executors.

## Class Diagram

```mermaid
--8<-- "class-diagrams/nextflow.executor.mmd"
```

!!! note
    Some classes may be excluded from the above diagram for brevity.

## Notes

The `Executor` class is the base class for all Nextflow executors. The main purpose of an `Executor` is to submit tasks to a scheduler. It uses a `TaskMonitor` to manage the lifecycle of all tasks and `TaskHandler`s to manage each individual task. See [nextflow.processor](nextflow.processor.md) for more details about these classes.

The built-in executors include the local executor (`LocalExecutor`) and the various grid executors (SLURM, PBS, LSF, etc), all of which extend `AbstractGridExecutor`. The `LocalExecutor` implements both "local" tasks (processes with a `script` or `shell` block) and "native" tasks (processes with an `exec` block).

The `BashWrapperBuilder` is used by executors to generate the wrapper script (`.command.run`) for a task, from a template script called `command-run.txt`, as well as the task configuration and the execution environment.
