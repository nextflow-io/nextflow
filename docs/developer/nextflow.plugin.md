
# `nextflow.plugin`

The `nextflow.plugin` package implements the plugin manager.

## Class Diagram

```{mermaid} diagrams/nextflow.plugin.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The plugin system uses the [PF4J](https://pf4j.org/) library, which allows for extension classes to be loaded at runtime. Each plugin includes a manifest of extension classes, all of which extend or implement some base class in Nextflow. The `Plugins` class can be used to query the available extensions for a given base class. Extensions can be assigned a priority using the `@Priority` annotation, to ensure that certain extensions are used over others when available.
