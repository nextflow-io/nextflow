
# `nextflow.extension`

The `nextflow.extension` package implements the channel operators.

## Class Diagram

```{mermaid} diagrams/nextflow.extension.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

Operators are implemented using the [GPars](http://gpars.org/1.2.1/guide/guide/dataflow.html) dataflow library. In general, an operator consumes one or more `DataflowReadChannel`s and produces one or more `DataflowWriteChannel`s.

See {ref}`operator-page` for details about each operator.
