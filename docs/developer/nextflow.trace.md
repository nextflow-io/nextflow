
# `nextflow.trace`

The `nextflow.trace` package defines the trace observer interface and implements several built-in trace observers.

## Class Diagram

```{mermaid} diagrams/nextflow.trace.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The `TraceObserver` interface defines a set of hooks into the workflow execution, such as when a workflow starts and completes, when a task starts and completes, and when an output file is published. The `Session` maintains a list of all observers and triggers each hook when the corresponding event occurs. Implementing classes can use these hooks to perform custom behaviors. In fact, this interface is used to implemented several core features, including the various execution reports, DAG renderer, and the integration with Nextflow Tower.
