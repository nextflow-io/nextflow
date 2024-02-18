
# `nextflow.container`

The `nextflow.container` package implements the integration with container runtimes.

## Class Diagram

```{mermaid} diagrams/nextflow.container.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The `ContainerBuilder` class is the base class for all container runtimes supported by Nextflow. It produces the container wrapper command for a given task run.

Executors that support containerized tasks insert this wrapper command into the task wrapper script (`.command.run`). Executors that are *container-native*, i.e. that launch the task wrapper itself inside a container, don't need to generate a container wrapper command.
