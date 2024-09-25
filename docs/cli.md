(cli-page)=

# Command line

Nextflow provides a robust command line interface (CLI) for the management and execution pipelines.

Simply run `nextflow` with no options or `nextflow -h` to see the list of available top-level options and commands. See {ref}`cli-reference` for the full list of subcommands with examples.

:::{note}
Nextflow options use a single dash prefix, e.g. `-foo`. Do not confuse with double dash notation, e.g. `--foo`, which is instead used for {ref}`Pipeline parameters <cli-params>`.
:::

## Basic usage

### Hard configuration override

Use the specified configuration file(s) overriding any defaults.

```console
$ nextflow -C my.config COMMAND [arg...]
```

The `-C` option is used to override *all* settings specified in the default config file. For soft override, please refer the `-c` option.

- Override **any** default configuration with a custom configuration file:

  ```console
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

  ```console
  $ nextflow -v
  nextflow version 20.07.1.5412
  ```

- The full version info with citation and website link:

  ```console
  $ nextflow -version
  N E X T F L O W
  version 20.07.1 build 5412
  created 24-07-2020 15:18 UTC (20:48 IDT)
  cite doi:10.1038/nbt.3820
  http://nextflow.io
  ```

## Running pipelines

The main purpose of the Nextflow CLI is to run Nextflow pipelines with the `run` command. Nextflow can execute a local script (e.g. `./main.nf`) or a remote project (e.g. `github.com/foo/bar`).

### Launching a remote project

To launch the execution of a pipeline project, hosted in a remote code repository, you simply need to specify its qualified name or the repository URL after the `run` command. The qualified name is formed by two parts: the `owner` name and the `repository` name separated by a `/` character.

In other words if a Nextflow project is hosted, for example, in a GitHub repository at the address `http://github.com/foo/bar`, it can be executed by entering the following command in your shell terminal:

```bash
nextflow run foo/bar
```

or using the project URL:

```bash
nextflow run http://github.com/foo/bar
```

If the project is found, it will be automatically downloaded to the Nextflow home directory (`$HOME/.nextflow` by default) and cached for subsequent runs.

:::{note}
You must use the `-hub` option to specify the hosting service if your project is hosted on a service other than GitHub, e.g. `-hub bitbucket`. However, the `-hub` option is not required if you use the project URL.
:::

Try this feature by running the following command:

```bash
nextflow run nextflow-io/hello
```

It will download a trivial example from the repository published at <http://github.com/nextflow-io/hello> and execute it on your computer.

If the `owner` is omitted, Nextflow will search your cached pipelines for a pipeline that matches the name specified. If no pipeline is found, Nextflow will try to download it using the `organization` name defined by the `NXF_ORG` environment variable (`nextflow-io` by default ).

:::{tip}
To access a private repository, specify the access credentials using the `-user` command line option. Then follow the interactive prompts to enter your password. Alternatively, define your private repository access credentials using Git. See {ref}`Git configuration <git-page>` for more information.
:::

### Using a specific revision

Any Git branch, tag, or commit of a project repository can be used when launching a pipeline by specifying the `-r` option:

```console
$ nextflow run nextflow-io/hello -r mybranch

or

```console
$ nextflow run nextflow-io/hello -r v1.1
```

These commands will execute two different project revisions based on the given Git branch/tag/commit.

(cli-params)=

### Pipeline parameters

Pipeline scripts can use an arbitrary number of parameters that can be overridden using the command line or Nextflow configuration files. Any script parameter can be specified on the command line by prefixing the parameter name with double-dash characters. For example:

```console
$ nextflow run <pipeline> --foo Hello
```

Then, the parameter can be accessed in the pipeline script using the `params.foo` identifier.

:::{note}
When the parameter name is formatted using `camelCase`, a second parameter is created with the same value using `kebab-case`, and vice versa.
:::

:::{warning}
When a command line parameter includes one or more glob characters, i.e. wildcards like `*` or `?`, the parameter value must be enclosed in quotes to prevent Bash expansion and preserve the glob characters. For example:

```console
$ nextflow run <pipeline> --files "*.fasta"
```
:::

## Managing projects

Nextflow seamlessly integrates with popular Git providers, including [BitBucket](http://bitbucket.org/), [GitHub](http://github.com), and [GitLab](http://gitlab.com) for managing Nextflow pipelines as version-controlled Git repositories.

The following commands allow you to perform basic operations to manage your projects.

### Listing available projects

The `list` command allows you to list all the projects you have downloaded in your computer. For example:

```bash
nextflow list
```

This prints a list similar to the following:

```
cbcrg/ampa-nf
cbcrg/piper-nf
nextflow-io/hello
nextflow-io/examples
```

### Showing project information

By using the `info` command you can show information from a downloaded project. For example:

```console
$ nextflow info hello
project name: nextflow-io/hello
repository  : http://github.com/nextflow-io/hello
local path  : $HOME/.nextflow/assets/nextflow-io/hello
main script : main.nf
revisions   :
* master (default)
  mybranch
  v1.1 [t]
  v1.2 [t]
```

Starting from the top it shows: the project name; the Git repository URL; the local directory where the project has been downloaded; the script that is executed when launched; the list of available revisions i.e. branches and tags. Tags are marked with a `[t]` on the right and the checked-out revision is marked with a `*` on the left.

### Pulling or updating a project

The `pull` command allows you to download a project from a GitHub repository or to update it if that repository has already been downloaded. For example:

```console
$ nextflow pull nextflow-io/hello
```

Alternatively, you can use the repository URL as the name of the project to pull:

```console
$ nextflow pull https://github.com/nextflow-io/hello
```

Downloaded pipeline projects are stored in your directory `$HOME/.nextflow/assets` directory.

### Viewing the project code

The `view` command shows the content of the pipeline script you have pulled. For example:

```console
$ nextflow view nextflow-io/hello
```

By adding the `-l` option to the example above it will list the content of the repository.

### Cloning a project into a directory

The `clone` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example:

```console
$ nextflow clone nextflow-io/hello target-dir
```

If the destination directory is omitted the specified project is cloned to a directory with the same name as the pipeline base name (e.g. `hello`) in the current directory.

The `clone` command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push back your changes by using the usual Git/GitHub workflow.

### Deleting a project

Downloaded pipelines can be deleted by using the `drop` command, as shown below:

```bash
nextflow drop nextflow-io/hello
```
