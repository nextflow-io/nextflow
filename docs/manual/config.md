(config-page)=

# Configuration

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations. Since each configuration file can contain conflicting settings, the sources are ranked to determine which settings are applied. Possible configuration sources, in order of priority:

1. Parameters specified on the command line (`--something value`)
2. Parameters provided using the `-params-file` option
3. Config file specified using the `-c my_config` option
4. The config file named `nextflow.config` in the current directory
5. The config file named `nextflow.config` in the workflow project directory
6. The config file `$HOME/.nextflow/config`
7. Values defined within the pipeline script itself (e.g. `main.nf`)

When more than one of these options for specifying configurations are used, they are merged, so that the settings in the first override the same settings appearing in the second, and so on.

:::{tip}
If you want to ignore any default configuration files and use only a custom one, use `-C <config file>`.
:::

## Syntax

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax:

```groovy
name = value
```

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (`true`, `false`) do not. Also note that values are typed. This means that, for example, `1` is different from `'1'` — the former is interpreted as the number one, while the latter is interpreted as a string value.

### Variables

Configuration properties can be used as variables in the configuration file by using the usual `$propertyName` or `${expression}` syntax.

For example:

```groovy
propertyOne = 'world'
anotherProp = "Hello $propertyOne"
customPath = "$PATH:/my/app/folder"
```

Please note, the usual rules for {ref}`string-interpolation` are applied, thus a string containing a variable reference must be wrapped in double-quote chars instead of single-quote chars.

The same mechanism allows you to access environment variables defined in the hosting system. Any variable name not defined in the Nextflow configuration file(s) is interpreted to be a reference to an environment variable with that name. So, in the above example, the property `customPath` is defined as the current system `PATH` to which the string `/my/app/folder` is appended.

### Comments

Configuration files use the same conventions for comments used by the Groovy or Java programming languages. Thus, use `//` to comment a single line, or `/*` .. `*/` to comment a block on multiple lines.

### Includes

A configuration file can include one or more configuration files using the keyword `includeConfig`. For example:

```groovy
process.executor = 'sge'
process.queue = 'long'
process.memory = '10G'

includeConfig 'path/foo.config'
```

When a relative path is used, it is resolved against the actual location of the including file.

## Config scopes

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope identifier, or grouping the properties in the same scope using the curly brackets notation. For example:

```groovy
alpha.x = 1
alpha.y = 'string value..'

beta {
     p = 2
     q = 'another string ..'
}
```

See {ref}`config-scopes` for the full list of config scopes.

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

(config-process)=

## Process configuration

The `process` scope allows you to specify default {ref}`directives <process-directives>` for processes in your pipeline.

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

:::{note}
The `withName` selector applies to a process even when it is included from a module under an alias. For example, `withName: hello` will apply to any process originally defined as `hello`, regardless of whether it is included under an alias. Similarly, it will not apply to any process not originally defined as `hello`, even if it is included under the alias `hello`.
:::

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

When mixing generic process configuration and selectors the following priority rules are applied (from lower to higher):

1. Process generic configuration.
2. Process specific directive defined in the workflow script.
3. `withLabel` selector definition.
4. `withName` selector definition.

For example:

```groovy
process {
    cpus = 4
    withLabel: foo { cpus = 8 }
    withName: bar { cpus = 32 }
}
```

Using the above configuration snippet, all workflow processes use 4 cpus if not otherwise specified in the workflow script. Moreover processes annotated with the `foo` label use 8 cpus. Finally the process named `bar` uses 32 cpus.

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
