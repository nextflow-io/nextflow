
# `nextflow.dag`

The `nextflow.dag` package implements the workflow DAG and renderers for several diagram formats.

## Class Diagram

```{mermaid} diagrams/nextflow.dag.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The workflow DAG defines the network of processes, channels, and operators that comprise a workflow. It is produced by the execution of the Nextflow script. See [nextflow.script](nextflow.script.md) for more details.

Implementations of the `DagRenderer` interface define how to render the workflow DAG to a particular diagram format. See {ref}`dag-visualisation` for more details.
