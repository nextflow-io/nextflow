
# `nextflow`

The `nextflow` package contains various top-level classes.

## Class Diagram

```{mermaid} diagrams/nextflow.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The `Nextflow` class implements several methods that are exposed to Nextflow scripts. See {ref}`implicit-functions` for details.

The `Channel` class implements the channel factory methods, and it is exposed directly to Nextflow scripts. See {ref}`channel-page` for details.

The `Session` class is the top-level representation of a Nextflow run, or "session". See [nextflow.script](nextflow.script.md) for more details about how a `Session` is created.
