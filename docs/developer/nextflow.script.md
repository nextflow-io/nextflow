
# `nextflow.script`

The `nextflow.script` package implements the parsing and execution of Nextflow scripts.

## Class Diagram

```{mermaid} diagrams/nextflow.script.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The execution of a Nextflow pipeline occurs in two phases. In the first phase, Nextflow parses and runs the script (using the language extensions in [nextflow.ast](nextflow.ast.md) and [nextflow.extension](nextflow.extension.md)), which constructs the workflow DAG. In the second phase, Nextflow executes the workflow.

```{note}
In DSL1, there was no separation between workflow construction and execution -- dataflow operators were executed as soon as they were constructed. DSL2 introduced lazy execution in order to separate process definition from execution, and thereby facilitate subworkflows and modules.
```
