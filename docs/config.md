(config-page)=

# Configuration

## Configuration file

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations. Since each configuration file may contain conflicting settings, they are applied in the following order (from lowest to highest priority):

1. Parameters defined in pipeline scripts (e.g. `main.nf`)
2. The config file `$HOME/.nextflow/config`
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

`NXF_ANSI_SUMMARY`
: Enables/disables ANSI completion summary: `true\|false` (default: print summary if execution last more than 1 minute).

`NXF_ASSETS`
: Defines the directory where downloaded pipeline repositories are stored (default: `$NXF_HOME/assets`)

`NXF_CACHE_DIR`
: :::{versionadded} 24.02.0-edge
  :::
: Defines the base cache directory when using the default cache store (default: `"$launchDir/.nextflow"`).

`NXF_CHARLIECLOUD_CACHEDIR`
: Directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CLASSPATH`
: Allows the extension of the Java runtime classpath with extra JAR files or class folders.

`NXF_CLOUDCACHE_PATH`
: :::{versionadded} 23.07.0-edge
  :::
: Defines the base cache path when using the cloud cache store.

`NXF_CLOUD_DRIVER`
: Defines the default cloud driver to be used if not specified in the config file or as command line option, either `aws` or `google`.

`NXF_CONDA_CACHEDIR`
: Directory where Conda environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CONDA_ENABLED`
: :::{versionadded} 22.08.0-edge
  :::
: Enable the use of Conda recipes defined by using the {ref}`process-conda` directive. (default: `false`).

`NXF_DEFAULT_DSL`
: :::{versionadded} 22.03.0-edge
  :::
: Defines the DSL version that should be used in not specified otherwise in the script of config file (default: `2`)

`NXF_DISABLE_CHECK_LATEST`
: :::{versionadded} 23.09.0-edge
  :::
: Nextflow automatically checks for a newer version of itself unless this option is enabled (default: `false`).

`NXF_DISABLE_JOBS_CANCELLATION`
: :::{versionadded} 21.12.0-edge
  :::
: Disables the cancellation of child jobs on workflow execution termination.

`NXF_DISABLE_PARAMS_TYPE_DETECTION`
: :::{versionadded} 23.07.0-edge
  :::
: Disables the automatic type detection of command line parameters.

`NXF_DISABLE_WAVE_SERVICE`
: :::{versionadded} 23.08.0-edge
  :::
: Disables the requirement for Wave service when enabling the Fusion file system.

`NXF_ENABLE_AWS_SES`
: :::{versionadded} 23.06.0-edge
  :::
: Enable to use of AWS SES native API for sending emails in place of legacy SMTP settings (default: `false`)

`NXF_ENABLE_FS_SYNC`
: :::{versionadded} 23.10.0
  :::
: When enabled the job script will execute Linux `sync` command on job completion. This may be useful to synchronize the job state over shared file systems (default: `false`)

`NXF_ENABLE_SECRETS`
: :::{versionadded} 21.09.0-edge
  :::
: Enable Nextflow secrets features (default: `true`)

`NXF_ENABLE_STRICT`
: :::{versionadded} 22.05.0-edge
  :::
: Enable Nextflow *strict* execution mode (default: `false`)

`NXF_ENABLE_VIRTUAL_THREADS`
: :::{versionadded} 23.05.0-edge
  :::
: :::{versionchanged} 23.10.0
  Enabled by default when using Java 21 or later.
  :::
: Enable the use of virtual threads in the Nextflow runtime (default: `false`)

`NXF_EXECUTOR`
: Defines the default process executor, e.g. `sge`

`NXF_FILE_ROOT`
: :::{versionadded} 23.05.0-edge
  :::
: The file storage path against which relative file paths are resolved.
: For example, with `NXF_FILE_ROOT=/some/root/path`, the use of `file('foo')` will be resolved to the absolute path `/some/root/path/foo`. A remote root path can be specified using the usual protocol prefix, e.g. `NXF_FILE_ROOT=s3://my-bucket/data`. Files defined using an absolute path are not affected by this setting.

`NXF_HOME`
: Nextflow home directory (default: `$HOME/.nextflow`).

`NXF_JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the `JAVA_HOME` variable if defined.

`NXF_JVM_ARGS`
: :::{versionadded} 21.12.1-edge
  :::
: Allows the setting Java VM options. This is similar to `NXF_OPTS` however it's only applied the JVM running Nextflow and not to any java pre-launching commands.

`NXF_LOG_FILE`
: The filename of the Nextflow log (default: `.nextflow.log`)

`NXF_OFFLINE`
: When `true` prevents Nextflow from automatically downloading and updating remote project repositories (default: `false`).
: :::{versionchanged} 23.09.0-edge
  This option also disables the automatic version check (see `NXF_DISABLE_CHECK_LATEST`).
  :::
: :::{versionchanged} 23.11.0-edge
  This option also prevents plugins from being downloaded. Plugin versions must be specified in offline mode, or else Nextflow will fail.
  :::

`NXF_OPTS`
: Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of `-Dkey[=value]` properties.

`NXF_ORG`
: Default `organization` prefix when looking for a hosted repository (default: `nextflow-io`).

`NXF_PARAMS_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the pipeline parameters file .

`NXF_PID_FILE`
: Name of the file where the process PID is saved when Nextflow is launched in background.

`NXF_PLUGINS_DEFAULT`
: Whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: `true`).

`NXF_PLUGINS_DIR`
: The path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins`).

`NXF_PLUGINS_TEST_REPOSITORY`
: :::{versionadded} 23.04.0
  :::
: Defines a custom plugin registry or plugin release URL for testing plugins outside of the main registry. See {ref}`testing-plugins` for more information.

`NXF_SCM_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the SCM config file .

`NXF_SINGULARITY_CACHEDIR`
: Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SINGULARITY_LIBRARYDIR`
: :::{versionadded} 21.09.0-edge
  :::
: Directory where remote Singularity images are retrieved. It should be a directory accessible to all compute nodes.

`NXF_SPACK_CACHEDIR`
: Directory where Spack environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SPACK_ENABLED`
: :::{versionadded} 23.02.0-edge
  :::
: Enable the use of Spack recipes defined by using the {ref}`process-spack` directive. (default: `false`).

`NXF_TEMP`
: Directory where temporary files are stored

`NXF_TRACE`
: Enable trace level logging for the specified packages. Equivalent to the `-trace` command-line option.

`NXF_VER`
: Defines which version of Nextflow to use.

`NXF_WORK`
: Directory where working files are stored (usually your *scratch* directory)

`NXF_WRAPPER_STAGE_FILE_THRESHOLD`
: :::{versionadded} 23.05.0-edge
  :::
: Defines the minimum size of the `.command.run` staging script for it to be written to a separate `.command.stage` file (default: `'1 MB'`).
: This setting is useful for executors that impose a size limit on job scripts.

`JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow.

`JAVA_CMD`
: Defines the path location of the Java binary command used to launch Nextflow.

`HTTP_PROXY`
: Defines the HTTP proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `http://user:password@proxy-host.com:port`.
  :::

`HTTPS_PROXY`
: Defines the HTTPS proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `https://user:password@proxy-host.com:port`.
  :::

`FTP_PROXY`
: :::{versionadded} 21.06.0-edge
  :::
: Defines the FTP proxy server. Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `ftp://user:password@proxy-host.com:port`.

`NO_PROXY`
: Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.

(config-feature-flags)=

## Feature flags

Some features can be enabled using the `nextflow.enable` and `nextflow.preview` flags. These flags can be specified in the pipeline script or the configuration file, and they are generally used to introduce experimental or other opt-in features.

`nextflow.enable.configProcessNamesValidation`

: When `true`, prints a warning for every `withName:` process selector that doesn't match a process in the pipeline (default: `true`).

`nextflow.enable.dsl`

: Defines the DSL version to use (`1` or `2`).

: :::{versionadded} 22.03.0-edge
  DSL2 is the default DSL version.
  :::

: :::{versionadded} 22.12.0-edge
  DSL1 is no longer supported.
  :::

`nextflow.enable.moduleBinaries`

: When `true`, enables the use of modules with binary scripts. See {ref}`module-binaries` for more information.

`nextflow.enable.strict`

: :::{versionadded} 22.05.0-edge
  :::

: When `true`, the pipeline is executed in "strict" mode, which introduces the following rules:

  - When reading a params file, Nextflow will fail if a dynamic param value references an undefined variable

  - When merging params from a config file with params from the command line, Nextflow will fail if a param is specified from both sources but with different types

  - When using the `join` operator, the `failOnDuplicate` option is `true` regardless of any user setting

  - When using the `join` operator, the `failOnMismatch` option is `true` (unless `remainder` is also `true`) regardless of any user setting

  - When using the `publishDir` process directive, the `failOnError` option is `true` regardless of any user setting

  - In a process definition, Nextflow will fail if an input or output tuple has only one element

  - In a process definition, Nextflow will fail if an output emit name is not a valid identifier (i.e. it should match the pattern `/[A-Za-z_][A-Za-z0-9_]*/`)

  - During a process execution, Nextflow will fail if a received input tuple does not have the same number of elements as was declared

  - During a process execution, Nextflow will fail if the `storeDir` directive is used with non-file outputs

  - Nextflow will fail if a pipeline param is referenced before it is defined

  - Nextflow will fail if multiple functions and/or processes with the same name are defined in a module script

`nextflow.preview.recursion`

: :::{versionadded} 21.11.0-edge
  :::

: *Experimental: may change in a future release.*

: When `true`, enables process and workflow recursion. See [this GitHub discussion](https://github.com/nextflow-io/nextflow/discussions/2521) for more information.

`nextflow.preview.topic`

: :::{versionadded} 23.11.0-edge
  :::

: *Experimental: may change in a future release.*

: When `true`, enables {ref}`topic channels <channel-topic>` feature.

## Strict config mode

:::{versionadded} 24.10.0
:::

Strict config mode is an alternative implementation of the config parser in Nextflow which enforces a stricter syntax. Whereas the existing parser interprets config files as Groovy scripts, which makes it overly powerful and error-prone, the strict config mode allows only a strict subset of Groovy syntax which is relevant to configuration.

Strict config mode can be enabled by setting `NXF_ENABLE_STRICT_CONFIG=true`.

The following statements are allowed in strict config mode:

- Assignment (e.g. `foo.bar = '3'`)
- Block (e.g. `foo { bar { ... } }`)
- Include (e.g. `includeConfig 'foo.config'`)

When assigning a config option, the right-hand side can be any Groovy expression, including closures which may include Groovy statements. However, only a subset of Groovy is supported within this context.

Notable syntax that is supported:

- numbers
- strings
- lists
- maps
- function calls
- closures
- variables (within closures)
- if-else statements (within closures)

Notable syntax that is not supported:

- function definitions
- class definitions
- try-catch statements
- lambdas (use closures instead)

Other syntax restrictions:

- The implicit use of environment variables is no longer supported:
  ```groovy
  env.MY_VAR = "$MY_VAR"                // incorrect
  env.MY_VAR = System.getenv('MY_VAR')  // correct
  ```

- The implicit `it` variable in closures is no longer supported:
  ```groovy
  foo.bar = ['foo', 'bar'].collect { "--$it" }.join(' ')        // incorrect
  foo.bar = ['foo', 'bar'].collect { it -> "--$it" }.join(' ')  // correct
  ```
