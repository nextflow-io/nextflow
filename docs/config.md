(config-page)=

# Configuration

## Configuration file

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations. Since each configuration file may contain conflicting settings, they are applied in the following order (from lowest to highest priority):

1. Parameters defined in pipeline scripts (e.g. `main.nf`)
2. The config file `$HOME/.nextflow/config`, or `$NXF_HOME/.nextflow/config` when `NXF_HOME` is set (see [`NXF` prefixed variables](reference/env-vars.html#nextflow-settings)).
3. The config file `nextflow.config` in the project directory
4. The config file `nextflow.config` in the launch directory
5. Config file specified using the `-c <config-file>` option
6. Parameters specified in a params file (`-params-file` option)
7. Parameters specified on the command line (`--something value`)

When more than one of these options for specifying configurations are used, they are merged, so that the settings in the first override the same settings appearing in the second, and so on.

:::{tip}
You can use the `-C <config-file>` option to use a single configuration file and ignore all other files.
:::

(config-syntax)=

## Syntax

The Nextflow configuration syntax is based on the Nextflow script syntax. It is designed for setting configuration options in a declarative manner while also allowing for dynamic expressions where appropriate.

A Nextflow config file may consist of any number of *assignments*, *blocks*, and *includes*. Config files may also contain comments in the same manner as scripts.

See {ref}`syntax-page` for more information about the Nextflow script syntax.

### Assignments

A config assignment consists of a config option and an expression separated by an equals sign:

```groovy
workDir = 'work'
docker.enabled = true
process.maxErrors = 10
```

A config option consists of an *option name* prefixed by any number of *scopes* separated by dots. Config scopes are used to group related config options. See {ref}`config-options` for the full set of config options.

The expression is typically a literal value such as a number, boolean, or string. However, any expression can be used:

```groovy
params.helper_file = "${projectDir}/assets/helper.txt"
```

### Blocks

A config scope can also be specified as a block, which may contain multiple configuration options. For example:

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

As a result, deeply nested config options can be assigned in various ways. For example, the following three assignments are equivalent:

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

A configuration file can include any number of other configuration files using the `includeConfig` keyword:

```groovy
process.executor = 'sge'
process.queue = 'long'
process.memory = '10G'

includeConfig 'path/foo.config'
```

Relative paths are resolved against the location of the including file.

:::{note}
Config includes can also be specified within config blocks. However, config files should only be included at the top level or in a [profile](#config-profiles) so that the included config file is valid on its own and in the context in which it is included.
:::

## Constants

The following constants are globally available in a Nextflow configuration file:

`baseDir`
: :::{deprecated} 20.04.0
  :::
: Alias for `projectDir`.

`launchDir`
: The directory where the workflow was launched.

`projectDir`
: The directory where the main script is located.

## Functions

The following functions are globally available in a Nextflow configuration file:

`env( name )`
: :::{versionadded} 24.11.0-edge
  :::
: Get the value of the environment variable with the specified name in the Nextflow launch environment.

(config-params)=

## Parameters

Pipeline parameters can be defined in the config file using the `params` scope:

```groovy
params.custom_param = 123
params.another_param = 'string value .. '

params {
    alpha_1 = true
    beta_2 = 'another string ..'
}
```

See {ref}`cli-params` for information about how to modify these on the command line.

(config-process)=

## Process configuration

The `process` scope allows you to specify {ref}`process directives <process-reference>` separately from the pipeline code.

For example:

```groovy
process {
    executor = 'sge'
    queue = 'long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
}
```

By using this configuration, all processes in your pipeline will be executed through the SGE cluster, with the specified settings.

(config-process-selectors)=

### Process selectors

The `withLabel` selectors allow the configuration of all processes annotated with a {ref}`process-label` directive as shown below:

```groovy
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
        queue = 'long'
    }
}
```

The above configuration example assigns 16 cpus, 64 Gb of memory and the `long` queue to all processes annotated with the `big_mem` label.

In the same manner, the `withName` selector allows the configuration of a specific process in your pipeline by its name. For example:

```groovy
process {
    withName: hello {
        cpus = 4
        memory = 8.GB
        queue = 'short'
    }
}
```

The `withName` selector applies both to processes defined with the same name and processes included under the same alias. For example, `withName: hello` will apply to any process originally defined as `hello`, as well as any process included under the alias `hello`.

Furthermore, selectors for the alias of an included process take priority over selectors for the original name of the process. For example, given a process defined as `foo` and included as `bar`, the selectors `withName: foo` and `withName: bar` will both be applied to the process, with the second selector taking priority over the first.

:::{tip}
Label and process names do not need to be enclosed with quotes, provided the name does not include special characters (`-`, `!`, etc) and is not a keyword or a built-in type identifier. When in doubt, you can enclose the label name or process name with single or double quotes.
:::

(config-selector-expressions)=

### Selector expressions

Both label and process name selectors allow the use of a regular expression in order to apply the same configuration to all processes matching the specified pattern condition. For example:

```groovy
process {
    withLabel: 'foo|bar' {
        cpus = 2
        memory = 4.GB
    }
}
```

The above configuration snippet sets 2 cpus and 4 GB of memory to the processes annotated with a label `foo` and `bar`.

A process selector can be negated prefixing it with the special character `!`. For example:

```groovy
process {
    withLabel: 'foo' { cpus = 2 }
    withLabel: '!foo' { cpus = 4 }
    withName: '!align.*' { queue = 'long' }
}
```

The above configuration snippet sets 2 cpus for the processes annotated with the `foo` label and 4 cpus to all processes *not* annotated with that label. Finally it sets the use of `long` queue to all process whose name does *not* start with `align`.

(config-selector-priority)=

### Selector priority

Process configuration settings are applied to a process in the following order (from lowest to highest priority):

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
    withLabel: foo { cpus = 8 }
    withName: bar { cpus = 16 }
    withName: 'baz:bar' { cpus = 32 }
}
```

With the above configuration:
- All processes will use 4 cpus (unless otherwise specified in their process definition).
- Processes annotated with the `foo` label will use 8 cpus.
- Any process named `bar` (or imported as `bar`) will use 16 cpus.
- Any process named `bar` (or imported as `bar`) invoked by a workflow named `baz` with use 32 cpus.

(config-profiles)=

## Config profiles

Configuration files can contain the definition of one or more *profiles*. A profile is a set of configuration attributes that can be selected during pipeline execution by using the `-profile` command line option.

Configuration profiles are defined by using the special scope `profiles`, which group the attributes that belong to the same profile using a common prefix. For example:

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

This configuration defines three different profiles: `standard`, `cluster`, and `cloud`, that each set different process
configuration strategies depending on the target runtime platform. The `standard` profile is used by default when no profile is specified.

:::{tip}
Multiple configuration profiles can be specified by separating the profile names with a comma, for example:

```bash
nextflow run <your script> -profile standard,cloud
```
:::

:::{danger}
When using the `profiles` feature in your config file, do NOT set attributes in the same scope both inside and outside a `profiles` context. For example:

```groovy
process.cpus = 1

profiles {
  foo {
    process.memory = '2 GB'
  }

  bar {
    process.memory = '4 GB'
  }
}
```

In the above example, the `process.cpus` attribute is not correctly applied because the `process` scope is also used in the `foo` and `bar` profiles.
:::

## Workflow handlers

Workflow event handlers can be defined in the config file, which is useful for handling pipeline events without having to modify the pipeline code:

```groovy
workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}

workflow.onError = {
    println "Error: something when wrong"
}
```

See {ref}`workflow-handlers` for more information.
