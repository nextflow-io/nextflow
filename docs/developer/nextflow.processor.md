
# `nextflow.processor`

The `nextflow.processor` package implements the execution and monitoring of tasks.

## Class Diagram

```{mermaid} diagrams/nextflow.processor.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

While the [`executor`](nextflow.executor.md) package defines how tasks are submitted to a particular execution environment (such as an HPC scheduler), the `processor` package defines how tasks are created and executed. As such, these packages work closely together, and in fact several components of the `Executor` interface, specifically the `TaskHandler` and `TaskMonitor`, are defined in this package.

The `TaskProcessor` is by far the largest and most complex class in this package. It implements both the dataflow operator for a given process as well as the task execution logic. In other words, it defines the mapping from an abstract process definition with concrete channel inputs into concrete task executions.

A `TaskRun` represents a particular task execution. There is also `TaskBean`, which is a serializable representation of a task. Legends say that `TaskBean` was originally created to support a "daemon" mode in which Nextflow would run on both the head node and the worker nodes, so the Nextflow "head" would need to send tasks to the Nextflow "workers". This daemon mode was never completed, but echoes of it remain (see `CmdNode`, `DaemonLauncher`, and the `nf-ignite` plugin).
