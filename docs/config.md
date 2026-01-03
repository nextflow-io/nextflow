(config-page)=

# Configuration

## Configuration file

When you launch a pipeline script, Nextflow looks for configuration files in multiple locations. Nextflow applies conflicting settings in the following order (from lowest to highest):

1. `$HOME/.nextflow/config`, or `$NXF_HOME/config` when {ref}`NXF_HOME <nxf-env-vars>` is set
2. `nextflow.config` in the project directory
3. `nextflow.config` in the launch directory
4. Config files specified with `-c <config-files>`

:::{tip}
You can use the `-C <config-file>` option to specify a fixed set of configuration files and ignore all other files.
:::

(config-syntax)=

## Syntax

Nextflow configuration uses the same syntax as Nextflow scripts. You can set configuration options declaratively and define dynamic expressions when needed.

Config files can contain any number of [assignments](#assignments), [blocks](#blocks), and [includes](#includes). You can add comments just like in scripts.

For more about information about the Nextflow script syntax, see {ref}`syntax-page`.

### Assignments

A config assignment sets a config option to an expression using an equals sign:

```groovy
workDir = 'work'
docker.enabled = true
process.maxErrors = 10
```

A config option has an *option name* with any number of *scopes* as prefixes, separated by dots. Scopes group related config options. For all options, see {ref}`config-options`.

Expressions are typically literal values (numbers, booleans, or strings). However, you can use any expression:

```groovy
params.helper_file = "${projectDir}/assets/helper.txt"
```

### Blocks

You can specify a config scope as a block with multiple options:

```groovy
// dot syntax
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'

// block syntax
docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}
```

You can assign deeply nested config options in multiple ways. The following three assignments are equivalent:

```groovy
executor.retry.maxAttempt = 5

executor {
    retry.maxAttempt = 5
}

executor {
    retry {
        maxAttempt = 5
    }
}
```

### Includes

You can include configuration files in other configuration files with the `includeConfig` keyword:

```groovy
process.executor = 'sge'
process.queue = 'long'
process.memory = '10G'

includeConfig 'path/extra.config'
```

Nextflow resolves relative paths from the including file's location.

:::{note}
You can specify config includes within config blocks. However, you should only include config files at the top level or in a [profile](#config-profiles). This ensures the included config file is valid on its own and in the context in which it is included.
:::

(config-constants)=

## Constants and functions

Nextflow configuration files have access to constants and functions from the global namespace. Constants allow you to reference runtime paths like project and launch directories, while functions enable dynamic operations such as reading environment variables.

For example, you can use `projectDir` constant to reference files relative to your project location:

```groovy
params.helper_file = "${projectDir}/assets/helper.txt"
```

Or, use the `env()` function to read environment variables:

```groovy
process.queue = env('MY_QUEUE')
```

For the full list of available constants and functions, see {ref}`stdlib-namespaces-global`.

(config-params)=

## Parameters

You can define pipeline parameters in config files using the `params` scope:

```groovy
// dot syntax
params.max_cpus = 64
params.publish_mode = 'copy'

// block syntax
params {
    max_cpus = 64
    publish_mode = 'copy'
}
```

:::{note}
When including a config file, Nextflow evaluates the included config with parameters defined before the include. Parameters defined after the include are not visible to the included config.
:::

You should declare parameters in the config file only when other config options use them. When a parameter is used in the script, you should declare it there and override it in config profiles as needed:

```nextflow
// main.nf
params.input = null
```

```groovy
// nextflow.config
params {
    publish_mode = 'copy'
}

workflow.output.mode = params.publish_mode

profiles {
    test {
        params.input = "${projectDir}/test/input.txt"
    }
}
```

See {ref}`cli-params` for information about how pipeline parameters are resolved at runtime.

(config-process)=

## Process configuration

You can use the `process` scope to specify {ref}`process directives <process-reference>` separately from the pipeline code.

For example:

```groovy
process {
    executor = 'sge'
    queue = 'long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
}
```

This configuration executes all processes through the SGE cluster with the specified settings.

(config-process-selectors)=

### Process selectors

You can use the `withLabel` selector to configure all processes annotated with a {ref}`process-label` directive:

```groovy
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
        queue = 'long'
    }
}
```

This configuration assigns 16 CPUs, 64 GB of memory, and the `long` queue to all processes with the `big_mem` label.

You can use the `withName` selector to configure a specific process by its name:

```groovy
process {
    withName: hello {
        cpus = 4
        memory = 8.GB
        queue = 'short'
    }
}
```

The `withName` selector matches both:

- Processes defined with that name
- Processes included under that alias

When a process is included with an alias, selectors for the alias take priority over selectors for the original name. For example, if you define a process as `hello` and include it as `sayHello`, both `withName: hello` and `withName: sayHello` apply, with `sayHello` taking priority.

:::{tip}
You don't need to enclose label and process names in quotes unless they contain special characters (`-`, `!`, etc.) or are keywords or built-in type identifiers. When in doubt, use single or double quotes.
:::

(config-selector-expressions)=

### Selector expressions

You can use regular expressions in label and process name selectors to apply the same configuration to all processes matching the pattern:

```groovy
process {
    withLabel: 'hello|bye' {
        cpus = 2
        memory = 4.GB
    }
}
```

This configuration requests 2 CPUs and 4 GB of memory for processes labeled as `hello` or `bye`.

You can negate a process selector by prefixing it with the special character `!`:

```groovy
process {
    withLabel: 'hello' { cpus = 2 }
    withLabel: '!hello' { cpus = 4 }
    withName: '!align.*' { queue = 'long' }
}
```

This configuration sets 2 CPUs for processes labeled as `hello` and 4 CPUs for all processes *not* labeled as `hello`. It also specifies the `long` queue for processes whose name does *not* start with `align`.

(config-selector-priority)=

### Selector priority

Nextflow applies process configuration settings in the following order (from lowest to highest priority):

1. Process configuration settings (without a selector)
2. Process directives in the process definition
3. `withLabel` selectors matching any of the process labels
4. `withName` selectors matching the process name
5. `withName` selectors matching the process included alias
6. `withName` selectors matching the process fully qualified name

For example:

```groovy
process {
    cpus = 4
    withLabel: hello { cpus = 8 }
    withName: bye { cpus = 16 }
    withName: 'aloha:bye' { cpus = 32 }
}
```

This configuration:

- Sets 4 CPUs for all processes (unless otherwise specified in their process definition)
- Sets 8 CPUs for processes annotated with the `hello` label
- Sets 16 CPUs for any process named `bye` (or imported as `bye`)
- Sets 32 CPUs for any process named `bye` (or imported as `bye`) invoked by a workflow named `aloha`

(config-profiles)=

## Config profiles

Configuration files can define one or more *profiles*. A profile is a set of configuration settings that can be selected during pipeline execution using the `-profile` command line option.

Configuration profiles are defined in the `profiles` scope. For example:

```groovy
profiles {
    standard {
        process.executor = 'local'
    }

    cluster {
        process.executor = 'sge'
        process.queue = 'long'
        process.memory = '10GB'
    }

    cloud {
        process.executor = 'cirrus'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }
}
```

This configuration defines three profiles: `standard`, `cluster`, and `cloud`. Each profile provides a different configuration for a given execution environment. Nextflow uses the `standard` profile by default when no profile is specified.

You can specify configuration profiles at runtime as a comma-separated list:

```bash
nextflow run main.nf -profile standard,cloud
```

Nextflow applies config profiles in the order in which they were defined in the config file, regardless of the order you specify them on the command line.

:::{versionadded} 25.02.0-edge
When using the {ref}`strict config syntax <updating-config-syntax>`, Nextflow applies profiles in the order you specify them on the command line.
:::

:::{danger}
When defining a profile in the config file, avoid using both the dot and block syntax for the same scope. For example:

```groovy
profiles {
    cluster {
        process.memory = '2 GB'
        process {
            cpus = 2
        }
    }
}
```

Due to a limitation of the legacy config parser, the second setting overwrites the first:

```console
$ nextflow config -profile cluster
process {
   cpus = 2
}
```

You can avoid this limitation by using the {ref}`strict config syntax <updating-config-syntax>`.
:::

(config-workflow-handlers)=

## Workflow handlers

:::{deprecated} 25.10.0
Use a {ref}`trace observer <plugins-trace-observers>` in a plugin to add custom workflow handlers to a pipeline via configuration.
:::

You can define workflow event handlers in the config file:

```groovy
workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}

workflow.onError = {
    println "Error: something went wrong"
}
```

This approach is useful for handling workflow events without modifying the pipeline code. See {ref}`workflow-handlers` for more information.
