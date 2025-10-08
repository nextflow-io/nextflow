# Plugin Spec

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2025-09-22
- Tags: plugins

## Summary 

Provide a way for external systems to understand key information about third-party plugins.

## Problem Statement

Nextflow plugins need a way to statically declare extensions to the Nextflow language so that external systems can extract information about a plugin without loading it in the JVM.

Primary use cases:

- The Nextflow language server needs to know about any config scopes, custom functions, etc, defined by a plugin, in order to recognize them in Nextflow scripts and config files.

- The Nextflow plugin registry (or other user interfaces) can use this information to provide API documentation.

## Goals or Decision Drivers

- External systems (e.g. language server) need to be able to understand plugins without having to load them in the JVM.

## Non-goals

- Defining specs for the core runtime and core plugins: these definitions are handled separately, although they may share some functionality with plugin specs.

## Considered Options

### Nextflow plugin system

Require external systems to use Nextflow's plugin system to load plugins at runtime in order to extract information about them.

- **Pro:** Allows any information to be extracted since the entire plugin is loaded

- **Con:** Requires the entire Nextflow plugin system to be reused or reimplemented. Not ideal for Java applications since the plugin system is implemented in Groovy, incompatible with non-JVM applications

- **Con:** Requires plugins to be downloaded, cached, loaded in the JVM, even though there is no need to use the plugin.

### Plugin spec

Define a plugin spec for every plugin release which is stored and served by the plugin registry as JSON.

- **Pro:** Allows any system to inspect plugin definitions through a standard JSON document, instead of downloading plugins and loading them into a JVM.

- **Con:** Requires the plugin spec to be generated at build-time and stored in the plugin registry.

- **Con:** Requires a standard format to ensure interoperability across different versions of Nextflow, the language server, and third-party plugins.

## Solution

Define a plugin spec for every plugin release which is stored and served by the plugin registry as JSON.

- Plugin developers only need to define [extension points](https://nextflow.io/docs/latest/plugins/developing-plugins.html#extension-points) as usual, and the Gradle plugin will extract the plugin spec and store it in the plugin registry as part of each plugin release.

- The language server can infer which third-party plugins are required from the `plugins` block in a config file. It will retrieve the appropriate plugin specs from the plugin registry.

A plugin spec consists of a list of *definitions*. Each definition has a *type* and a *spec*.

Example:

```json
{
    "$schema": "https://raw.githubusercontent.com/nextflow-io/schemas/main/plugin/v1/schema.json",
    "definitions": [
        {
            "type": "ConfigScope",
            "spec": {
                // ...
            }
        },
        {
            "type": "Function",
            "spec": {
                // ...
            }
        },
    ]
}
```

The following types of definitions are allowed:

**ConfigScope**

Defines a top-level config scope. The spec consists of a *name*, an optional *description*, and *children*.

The children should be a list of definitions corresponding to nested config scopes and options. The following definitions are allowed:

- **ConfigOption**: Defines a config option. The spec consists of a *description* and *type*.

- **ConfigScope**: Defines a nested config scope, using the same spec as for top-level scopes.

Example:

```json
{
    "type": "ConfigScope",
    "spec": {
        "name": "hello",
        "description": "The `hello` scope controls the behavior of the `nf-hello` plugin.",
        "children": [
            {
                "type": "ConfigOption",
                "spec": {
                    "name": "message",
                    "description": "Message to print to standard output when the plugin is enabled.",
                    "type": "String"
                }
            }
        ]
    }
}
```

**Factory**

Defines a channel factory that can be included in Nextflow scripts. The spec is the same as for functions.

**Function**

Defines a function that can be included in Nextflow scripts. The spec consists of a *name*, an optional *description*, a *return type*, and a list of *parameters*. Each parameter consists of a *name* and a *type*.

Example:

```json
{
    "type": "Function",
    "spec": {
        "name": "sayHello",
        "description": "Say hello to the given target",
        "returnType": "void",
        "parameters": [
            {
                "name": "target",
                "type": "String"
            }
        ]
    }
}
```

**Operator**

Defines a channel operator that can be included in Nextflow scripts. The spec is the same as for functions.

## Rationale & discussion

Now that there is a Gradle plugin for building Nextflow plugins and a registry to publish and retrieve plugins, it is possible to generate, publish, and retrieve plugin specs in a way that is transparent to plugin developers.

Plugins specs adhere to a pre-defined [schema](https://raw.githubusercontent.com/nextflow-io/schemas/main/plugin/v1/schema.json) to ensure consistency across different versions of Nextflow. In the future, new versions of the schema can be defined as needed to support new behaviors or requirements.
