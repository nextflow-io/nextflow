(cli-page)=

# Command line

Nextflow provides a robust command line interface for the management and execution pipelines.

Simply run `nextflow` with no options or `nextflow -h` to see the list of available top-level options and commands.

## Basic usage

### Hard configuration override

Use the specified configuration file(s) overriding any defaults.

```console
$ nextflow -C my.config COMMAND [arg...]
```

The `-C` option is used to override *all* settings specified in the default config file. For soft override, please refer the `-c` option.

- Override **any** default configuration with a custom configuration file:

  ```
  $ nextflow -C my.config run nextflow-io/hello
  ```

### JVM properties

Set JVM properties.

```console
$ nextflow -Dkey=value COMMAND [arg...]
```

This options allows the definition of custom Java system properties that can be used to properly configure or fine tuning the JVM instance used by the Nextflow runtime.

For specifying other JVM level options, please refer to the {ref}`config-env-vars` section.

- Add JVM properties to the invoked pipeline:

  ```console
  $ nextflow -Dfile.encoding=UTF-8 run nextflow-io/hello
  ```

### Execution as a background job

Execute `nextflow` in the background.

```console
$ nextflow -bg COMMAND [arg...]
```

The `-bg` option is used to invoke the nextflow execution in the background and allows the user to continue interacting with the terminal. This option is similar to `nohup` in behavior.

- Invoke any execution as a background job:

  ```console
  $ nextflow -bg run nextflow-io/hello
  ```

### Soft configuration override

Add the specified file to configuration set.

```console
$ nextflow -c nxf.config COMMAND [arg...]
```

The `-c` option is used to append a new configuration to the default configuration. The `-c` option allows us to update the config in an additive manner. For **hard override**, refer to the `-C` option.

- Update *some* fields of the default config for any pipeline:

  ```console
  $ nextflow -c nxf.config run nextflow-io/hello
  ```

### Docker driven execution

:::{deprecated} 23.09.0-edge
:::

Launch Nextflow via Docker.

```console
$ nextflow -dockerize COMMAND [arg...]
```

The `-dockerize` option is used to invoke the execution of Nextflow within a Docker container itself without installing a Java VM in the hosting environment.

This option is *not* needed to run containerised pipeline jobs. For invoking a pipeline with the `docker` profile or executor, please refer to the `-with-docker` options in the `run` command. When using the `-dockerize` option in combination with containerized tasks, Nextflow will launch the tasks as sibling containers in the host environment (i.e. no Docker-in-Docker).

- Invoke `nextflow` as a Docker container to execute a pipeline:

  ```console
  $ nextflow -dockerize run nextflow-io/hello
  ```

### Help

Print the help message.

```console
$ nextflow -h
```

The `-h` option prints out the overview of the CLI interface and enumerates the top-level *options* and *commands*.

### Execution logs

Sets the path of the nextflow log file.

```console
$ nextflow -log custom.log COMMAND [arg...]
```

The `-log` option takes a path of the new log file which to be used instead of the default `.nextflow.log` or to save logs files to another directory.

- Save all execution logs to the custom `/var/log/nextflow.log` file:

  ```console
  $ nextflow -log /var/log/nextflow.log run nextflow-io/hello
  ```

### Quiet execution

Disable the printing of information to the terminal.

```console
$ nextflow -q COMMAND [arg...]
```

The `-q` option suppresses the banner and process-related info, and exits once the execution is completed. Please note that it does not affect any explicit print statement within a pipeline.

- Invoke the pipeline execution without the banner and pipeline information:

  ```console
  $ nextflow -q run nextflow-io/hello
  ```

### Logging to a syslog server

Send logs to [Syslog](https://en.wikipedia.org/wiki/Syslog) server endpoint.

```console
$ nextflow -syslog localhost:1234 COMMAND [arg...]
```

The `-syslog` option is used to send logs to a Syslog logging server at the specified endpoint.

- Send the logs to a Syslog server at specific endpoint:

  ```console
  $ nextflow -syslog localhost:1234 run nextflow-io/hello
  ```

### Version

Print the Nextflow version information.

```console
$ nextflow -v
```

The `-v` option prints out information about Nextflow, such as the version and build. The `-version` option in addition prints out the citation reference and official website.

- The short version:

  ```
  $ nextflow -v
  nextflow version 20.07.1.5412
  ```

- The full version info with citation and website link:

  ```
  $ nextflow -version
  N E X T F L O W
  version 20.07.1 build 5412
  created 24-07-2020 15:18 UTC (20:48 IDT)
  cite doi:10.1038/nbt.3820
  http://nextflow.io
  ```

(cli-commands)=

## Commands

The Nextflow CLI provides subcommands for a variety of actions, including:

- `nextflow run`: run a pipeline
- `nextflow config`: view the resolved configuration of a pipeline
- `nextflow log`: view metadata of previous runs

Refer to the {ref}`cli-reference` for the full list of subcommands and examples for each.

(cli-params)=

## Pipeline parameters

Pipeline scripts can use an arbitrary number of parameters that can be overridden, either using the command line or the Nextflow configuration file. Any script parameter can be specified on the command line, prefixing the parameter name with double dash characters, e.g.:

```bash
nextflow run <my script> --foo Hello
```

Then, the parameter can be accessed in the pipeline script using the `params.foo` identifier.

:::{note}
When the parameter name is formatted using `camelCase`, a second parameter is created with the same value using `kebab-case`, and vice versa.
:::

:::{warning}
When a command line parameter includes one or more glob characters, i.e. wildcards like `*` or `?`, the parameter value must be enclosed in quotes to prevent Bash expansion and preserve the glob characters. For example:

```bash
nextflow run <my script> --files "*.fasta"
```
:::
