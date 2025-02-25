
# `nextflow.config`

The `nextflow.config` package contains the implementation of the Nextflow configuration.

## Class Diagram

```{mermaid} diagrams/nextflow.config.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

Any command that parses Nextflow config files (`config`, `run`, etc) uses the `ConfigBuilder` to build a `ConfigMap` from a set of config files. The `ConfigBuilder` itself uses a `ConfigParser` to parse the config files.

The Nextflow configuration language is essentially Groovy with some extensions. These extensions are implemented in `ConfigBase` and `ConfigTransformImpl`.
