(cli-page)=

# Command line interface (CLI)

Nextflow provides a robust command line interface for the management and execution pipelines.

Simply run `nextflow` with no options or `nextflow -h` to see the list of available top-level options and commands.

(cli-options)=

## Options

The top-level options are meant to be invoked in relation to the core Nextflow application and are applied to all commands. For options specific to any command, refer the CLI Commands section.

:::{note}
Nextflow options use a single dash prefix, e.g. `-foo`. Do not confuse with double dash notation, e.g. `--foo`, which is instead used for {ref}`Pipeline parameters <cli-params>`.
:::

Available options:

`-C`
: Use the specified configuration file(s) overriding any defaults.

`-D`
: Set JVM properties.

`-bg`
: Execute nextflow in background.

`-c, -config`
: Add the specified file to configuration set.

`-d, -dockerize`
: :::{deprecated} 23.09.0-edge
  :::
: Launch nextflow via Docker (experimental).

`-h`
: Print this help.

`-log`
: Set nextflow log file path.

`-q, -quiet`
: Do not print information messages.

`-remote-debug`
: Enable JVM interactive remote debugging (experimental).

`-syslog`
: Send logs to syslog server (e.g. localhost:514).

`-trace`
: Enable trace level logging for the specified packages. Multiple packages can be provided separating them with a comma, e.g. `-trace nextflow,io.seqera`.

`-v, -version`
: Print the program version.

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

(cli-clean)=

### clean

Clean up *cache* and *work* directories.

**Usage**

```console
$ nextflow clean [run_name|session_id] [options]
```

**Description**

Upon invocation within a directory, `nextflow` creates a project specific `.nextflow.log` file, `.nextflow` cache directory as well as a `work` directory. The `clean` command is designed to facilitate removal of these files from previous executions. A list of run names and session ids can be generated by invoking `nextflow log -q`.

If no run name or session id is provided, it will clean the latest run.

**Options**

`-after`
: Clean up runs executed *after* the specified one.

`-before`
: Clean up runs executed *before* the specified one.

`-but`
: Clean up all runs *except* the specified one.

`-n, -dry-run`
: Print names of files to be removed without deleting them.

`-f, -force`
: Force clean command.

`-h, -help`
: Print the command usage.

`-k, -keep-logs`
: Removes only temporary files but retains execution log entries and metadata.

`-q, -quiet`
: Do not print names of files removed.

**Examples**

Dry run to remove work directories for the run name `boring_euler`:

```console
$ nextflow clean boring_euler -n

Would remove work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
Would remove work/3f/70944c7a549b6221e1ccc7b4b21b62
Would remove work/0e/2ebdba85f76f6068b21a1bcbf10cab
```

Remove work directories for the run name `boring_euler`.

```console
$ nextflow clean boring_euler -f

Removed work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
Removed work/3f/70944c7a549b6221e1ccc7b4b21b62
Removed work/0e/2ebdba85f76f6068b21a1bcbf10cab
```

Remove the execution entries *except* for a specific execution.

```console
$ nextflow clean -but tiny_leavitt -f

Removed work/1f/f1ea9158fb23b53d5083953121d6b6
Removed work/bf/334115deec60929dc18edf0010032a
Removed work/a3/06521d75da296d4dd7f4f8caaddad8
```

Dry run to remove the execution data *before* a specific execution.

```console
$ nextflow clean -before tiny_leavitt -n

Would remove work/5d/ad76f7b7ab3500cf616814ef644b61
Would remove work/c4/69a82b080a477612ba8d8e4c27b579
Would remove work/be/a4fa2aa38f76fd324958c81c2e4603
Would remove work/54/39116773891c47a91e3c1733aad4de
```

Dry run to remove the execution data *after* a specific execution.

```console
$ nextflow clean -after focused_payne -n

Would remove work/1f/f1ea9158fb23b53d5083953121d6b6
Would remove work/bf/334115deec60929dc18edf0010032a
Would remove work/a3/06521d75da296d4dd7f4f8caaddad8
```

Dry run to remove the temporary execution data for a specific execution, while keeping the log files.

```console
$ nextflow clean -keep-logs tiny_leavitt -n

Would remove temp files from work/1f/f1ea9158fb23b53d5083953121d6b6
Would remove temp files from work/bf/334115deec60929dc18edf0010032a
Would remove temp files from work/a3/06521d75da296d4dd7f4f8caaddad8
```

### clone

Clone a remote project into a folder.

**Usage**

```console
$ nextflow clone [options] [project]
```

**Description**

The `clone` command downloads a pipeline from a Git-hosting platform into the *current directory* and modifies it accordingly. For downloading a pipeline into the global cache `~/.nextflow/assets`, please refer to the `nextflow pull` command.

**Options**

`-d, -deep`
: Create a shallow clone of the specified depth.

`-h, -help`
: Print the command usage.

`-hub` (`github`)
: Service hub where the project is hosted. Options: `gitlab` or `bitbucket`.

`-r` (`master`)
: Revision of the project to clone (either a git branch, tag or commit SHA number).

`-user`
: Private repository user name.

**Examples**

Clone the latest revision of a pipeline.

```console
$ nextflow clone nextflow-io/hello
nextflow-io/hello cloned to: hello
```

Clone a specific revision of a pipeline.

```console
$ nextflow clone nextflow-io/hello -r v1.1
nextflow-io/hello cloned to: hello
```

### config

Print the resolved pipeline configuration.

**Usage**

```console
$ nextflow config [options] [project name or path]
```

**Description**

The `config` command is used for printing the project's configuration i.e. the `nextflow.config` and is especially useful for understanding the resolved profiles and parameters that Nextflow will use run a pipeline. For in-depth information, please refer the {ref}`config-profiles` section.

**Options**

`-flat`
: Print config using flat notation.

`-h, -help`
: Print the command usage.

`-profile`
: Choose a configuration profile.

`-properties`
: Print config using Java properties notation.

`-r, -revision`
: Revision of the project (either a git branch, tag or commit SHA number).

`-a, -show-profiles`
: Show all configuration profiles.

`-sort`
: Sort config attributes.

`-value`
: :::{versionadded} 23.08.0-edge
  :::
: Print the value of a config option, or fail if the option is not defined.

**Examples**

Print out the inferred config using a the default group key-value notation.

```console
$ nextflow config

docker {
    enabled = true
}

process {
    executor = 'local'
}
```

Print out the config using a flat notation.

```console
$ nextflow config -flat

docker.enabled = true
process.executor = 'local'
```

Print out the config using the Java properties notation.

```console
$ nextflow config -properties

docker.enabled = true
process.executor = local
```

Print out the value of a specific configuration property.

```console
$ nextflow config -value process.executor
local
```

Print out all profiles from the project's configuration.

```console
$ nextflow config -show-profiles

docker {
    enabled = true
}

profiles {
    standard {
        process {
            executor = 'local'
        }
    }
    cloud {
        process {
            executor = 'cirrus'
            container = 'cbcrg/imagex'
        }
    }
}
```

(cli-console)=

### console

Launch the Nextflow interactive console.

**Usage**

```console
$ nextflow console
```

**Description**

The `console` command provides a Graphical User Interface (GUI) and an interactive REPL (Read-Eval-Print-Loop) for quick experimentation.

**Options**

None available

**Examples**

Launch the `console` GUI.

```console
$ nextflow console
```

### drop

Delete the local copy of a project.

**Usage**

```console
$ nextflow drop [options] [project]
```

**Description**

The `drop` command is used to remove the projects which have been downloaded into the global cache. Please refer the `list` command for generating a list of downloaded pipelines.

**Options**

`-a, -all-revisions`
: For specified project, drop all revisions.

`-f`
: Delete the repository without taking care of local changes.

`-h, -help`
: Print the command usage.

`-r, -revision`
: Revision of the project to drop (either a git branch, tag or commit SHA number).

**Examples**

Drop the `nextflow-io/hello` project.

```console
$ nextflow drop nextflow-io/hello
```

Forcefully drop the `nextflow-io/hello` pipeline, ignoring any local changes.

```console
$ nextflow drop nextflow-io/hello -f
```

### fs

Perform basic filesystem operations.

**Usage**

```console
$ nextflow fs [subcommands]
```

**Description**

The `fs` command is used to perform filesystem operations like copy, move, delete, list directory, etc. Like the `file()` method, it can work with local files, remote URLs, and remote object storage. Storage credentials can be provided through the same manner as launching a pipeline (Nextflow config, environment vars, etc).

**Options**

`-h, -help`
: Print the command usage.

**Examples**

List a directory.

```console
$ nextflow fs list <directory>
```

Print the contents of a file to standard output.

```console
$ nextflow fs cat <file>
```

Copy a file or directory.

```console
$ nextflow fs cp <source> <target>
```

Move a file or directory.

```console
$ nextflow fs mv <source> <target>
```

Delete a file or directory.

```console
$ nextflow fs rm <path>
```

:::{versionadded} 23.10.0
:::

Print file or directory attributes.

```console
$ nextflow fs stat <path>
```

### help

Print the top-level help or specific help for a command.

**Usage**

```console
$ nextflow help [options] [command]
```

**Description**

The `help` command prints out the overview of the CLI interface and enumerates the top-level *options* and *commands*. Note that this command is equivalent to simply invoking `nextflow` at the command line.

**Options**

`-h, -help`
: Print the command usage.

**Examples**

Invoke the `help` option for the `drop` command.

```console
$ nextflow help drop

Delete the local copy of a project
Usage: drop [options] name of the project to drop
   Options:
     -f
          Delete the repository without taking care of local changes
          Default: false
     -h, -help
          Print the command usage
          Default: false
```

### info

Print project or system runtime information.

**Usage**

```console
$ nextflow info [options] [project]
```

**Description**

The `info` command prints out the nextflow runtime information about the hardware as well as the software versions of the Nextflow version and build, operating system, and Groovy and Java runtime. It can also be used to display information about a specific project; in this case, note how revisions marked as `P` are pulled locally.

If no run name or session id is provided, it will clean the latest run.

**Options**

`-u, -check-updates`
: Check for remote updates.

`-d`
: Show detailed information.

`-h, -help`
: Print the command usage.

`-o` (`text`)
: Output format, either `text`, `json` or `yaml`.

**Examples**

Display Nextflow runtime and system info:

```console
$ nextflow info

  Version: 20.07.1 build 5412
  Created: 24-07-2020 15:18 UTC (20:48 IDT)
  System: Mac OS X 10.15.6
  Runtime: Groovy 2.5.11 on OpenJDK 64-Bit Server VM 1.8.0_192-b01
  Encoding: UTF-8 (UTF-8)
```

Display information about a specific project:

```console
$ nextflow info nextflow-io/hello

  project name: nextflow-io/hello
  repository  : https://github.com/nextflow-io/hello
  local path  : /Users/evanfloden/.nextflow/assets/nextflow-io/hello
  main script : main.nf
  revisions   :
  P master (default)
    mybranch
    testing
  P v1.1 [t]
    v1.2 [t]
```

### inspect

:::{versionadded} 23.09.0-edge
:::

Inspect process settings in a pipeline project. Currently only supports the `container` directive.

**Usage**

```console
$ nextflow inspect [options] [project]
```

**Description**

The `inspect` command allows you to determine the container for each process in a pipeline without running the pipeline. It prints to stdout a listing of containers for each process, formatted either as JSON or Nextflow configuration.

**Options**

`-concretize`
: Build the container images resolved by the inspect command.

`-format` (`json`)
: Inspect output format. Can be `json` or `config`.

`-i, -ignore-errors`
: Ignore errors while inspecting the pipeline.

`-params-file`
: Load script parameters from a JSON/YAML file.

`-profile`
: Use the given configuration profile(s).

`-r, revision`
: Revision of the project to inspect (either a git branch, tag or commit SHA number).

**Examples**

Get the list of containers used by a pipeline.

```console
$ nextflow inspect nextflow-io/hello
```

Specify parameters as with the `run` command:

```console
$ nextflow inspect main.nf --alpha 1 --beta foo
```

### kuberun

Launch a Nextflow pipeline on a Kubernetes cluster.

**Usage**

```console
$ nextflow kuberun [options] [project]
```

**Description**

The `kuberun` command builds upon the `run` command and offers a deep integration with the Kubernetes execution environment. This command deploys the Nextflow runtime as a Kubernetes pod and assumes that you've already installed the `kubectl` CLI. The `kuberun` command does not allow the execution of local Nextflow scripts. For more information please refer to the {ref}`k8s-page` page.

**Options**

The `kuberun` command supports the following options from [`run`](#run):

- `-cache`
- `-disable-jobs-cancellation`
- `-dsl1`
- `-dsl2`
- `-dump-channels`
- `-dump-hashes`
- `-e.<key>=<value>`
- `-entry`
- `-h, -help`
- `-hub`
- `-latest`
- `-main-script`
- `-name`
- `-offline`
- `-params-file`
- `-plugins`
- `-preview`
- `-process.<key>=<value>`
- `-profile`
- `-qs, -queue-size`
- `-resume`
- `-r, -revision`
- `-stub, -stub-run`
- `-user`
- `-with-conda`
- `-with-dag`
- `-N, -with-notification`
- `-with-report`
- `-with-spack`
- `-with-timeline`
- `-with-tower`
- `-with-trace`
- `-with-wave`
- `-with-weblog`
- `-without-spack`
- `-without-wave`
- `-w, -work-dir`

The following new options are also available:

`-head-cpus`
: :::{versionadded} 22.01.0-edge
  :::
: Specify number of CPUs requested for the Nextflow pod.

`-head-image`
: :::{versionadded} 22.07.1-edge
  :::
: Specify the container image for the Nextflow driver pod.

`-head-memory`
: :::{versionadded} 22.01.0-edge
  :::
: Specify amount of memory requested for the Nextflow pod.

`-head-prescript`
: :::{versionadded} 22.05.0-edge
  :::
: Specify script to be run before the Nextflow pod starts.

`-n, -namespace`
: Specify the K8s namespace to use.

`-remoteConfig`
: Add the specified file from the K8s cluster to configuration set.

`-remoteProfile`
: Choose a configuration profile in the remoteConfig.

`-v, -volume-mount`
: Volume claim mounts, e.g. `my-pvc:/mnt/path`.

**Examples**

Execute a pipeline into a Kubernetes cluster.

```console
$ nextflow kuberun nextflow-io/hello
```

### list

List all downloaded projects.

**Usage**

```console
$ nextflow list [options]
```

**Description**

The `list` commands prints a list of the projects which are already downloaded into the global cache `~/.nextflow/assets`.

**Options**

`-a, -all-revisions`
: For each project, also list revisions.

`-d`
: Show commit information for revisions (in conjunction with `-a`)).

`-h, -help`
: Print the command usage.

**Examples**

List the downloaded pipelines.

```console
$ nextflow list

nextflow-io/hello
nextflow-hub/fastqc
```

(cli-log)=

### log

Print the execution history and log information.

**Usage**

```console
$ nextflow log [options] [run_name | session_id]
```

**Description**

The `log` command is used to query the execution metadata associated with pipelines executed by Nextflow. The list of executed pipelines can be generated by running `nextflow log`. Instead of run name, it's also possible to use a session id. Moreover, this command contains multiple options to facilitate the queries and is especially useful while debugging a pipeline and while inspecting pipeline execution metadata.

**Options**

`-after`
: Show log entries for runs executed *after* the specified one.

`-before`
: Show log entries for runs executed *before* the specified one.

`-but`
: Show log entries for runs executed *but* the specified one.

`-f, -fields`
: Comma-separated list of fields to include in the printed log. Use the `-l` option to see the list of available fields.

`-F, -filter`
: Filter log entries by a custom expression, e.g. `process =~ /foo.*/ && status == 'COMPLETED'`.

`-h, -help`
: Print the command usage.

`-l, -list-fields`
: Show all available fields.

`-quiet`
: Show only run names.

`-s`
: Character used to separate column values.

`-t, -template`
: Text template used to each record in the log.

**Examples**

Listing the execution logs of previous invocations of all pipelines in a project.

```console
$ nextflow log

TIMESTAMP           DURATION        RUN NAME        STATUS  REVISION ID     SESSION ID                              COMMAND
2020-10-07 11:52:24 2.1s            focused_payne   OK      96eb04d6a4      af6adaaa-ad4f-48a2-9f6a-b121e789adf5    nextflow run nextflow-io/hello -r master
2020-10-07 11:53:00 3.1s            tiny_leavitt    OK      e3b475a61b      4d3b95c5-4385-42b6-b430-c865a70d56a4    nextflow run ./tutorial.nf
2020-10-07 11:53:29 2.5s            boring_euler    OK      e3b475a61b      a6276975-7173-4208-ae09-ab9d6dce8737    nextflow run tutorial.nf
```

Listing only the *run names* of the execution logs of all pipelines invocations in a project.

```console
$ nextflow log -quiet

focused_payne
tiny_leavitt
boring_euler
```

List the execution entries *only* a specific execution.

```console
$ nextflow log tiny_leavitt

work/1f/f1ea9158fb23b53d5083953121d6b6
work/bf/334115deec60929dc18edf0010032a
work/a3/06521d75da296d4dd7f4f8caaddad8
```

List the execution entries *after* a specific execution.

```console
$ nextflow log -after tiny_leavitt

work/92/c1a9cd9a96e0531d81ca69f5dc3bb7
work/3f/70944c7a549b6221e1ccc7b4b21b62
work/0e/2ebdba85f76f6068b21a1bcbf10cab
```

List the execution entries *before* a specific execution.

```console
$ nextflow log -before tiny_leavitt

work/5d/ad76f7b7ab3500cf616814ef644b61
work/c4/69a82b080a477612ba8d8e4c27b579
work/be/a4fa2aa38f76fd324958c81c2e4603
work/54/39116773891c47a91e3c1733aad4de
```

List the execution entries *except* for a specific execution.

```console
$ nextflow log -but tiny_leavitt

work/5d/ad76f7b7ab3500cf616814ef644b61
work/c4/69a82b080a477612ba8d8e4c27b579
work/be/a4fa2aa38f76fd324958c81c2e4603
work/54/39116773891c47a91e3c1733aad4de
```

Filter specific fields from the execution log of a process.

```console
$ nextflow log tiny_leavitt -f 'process,exit,hash,duration'

splitLetters        0       1f/f1ea91       112ms
convertToUpper      0       bf/334115       144ms
convertToUpper      0       a3/06521d       139ms
```

Filter fields from the execution log of a process based on a criteria.

```console
$ nextflow log tiny_leavitt -F 'process =~ /splitLetters/'

work/1f/f1ea9158fb23b53d5083953121d6b6
```

(cli-plugin)=

### plugin

Manage plugins and run plugin-specific commands.

```console
$ nextflow plugin <subcommand> [options]
```

The `plugin` command provides several subcommands for managing and using plugins:

`install <plugin[@version],..>`

: Install a plugin. Multiple plugins can be specified as a comma-separated list. Each plugin id consists of a name and optional version separated by a `@`.

`<plugin>:<subcommand> [options]`

: Execute a plugin-specific command.

### pull

Download or update a project.

**Usage**

```console
$ nextflow pull [options] [project]
```

**Description**

The `pull` command downloads a pipeline from a Git-hosting platform into the global cache `~/.nextflow/assets` and modifies it accordingly. For downloading a pipeline into a local directory, please refer to the `nextflow clone` command.

**Options**

`-a, -all`
: Update all downloaded projects.

`-d, -deep`
: Create a shallow clone of the specified depth.

`-h, -help`
: Print the command usage.

`-hub` (`github`)
: Service hub where the project is hosted. Options: `gitlab` or `bitbucket`

`-r, -revision`
: Revision of the project to pull (either a git branch, tag or commit SHA number).
: When passing a git tag or branch, the `workflow.revision` and `workflow.commitId` fields are populated. When passing only the commit hash, `workflow.revision` is not defined.

`-user`
: Private repository user name.

**Examples**

Download a new pipeline or pull the latest revision for a specific project.

```console
$ nextflow pull nextflow-io/hello

Checking nextflow-io/hello ...
done - revision: 96eb04d6a4 [master]
```

Pull the latest revision for all downloaded projects.

```console
$ nextflow pull -all

Checking nextflow-io/hello ...
done - revision: 96eb04d6a4 [master]
Checking nextflow-hub/fastqc ...
done - revision: 087659b18e [master]
```

Download a specific revision of a new project or pull the latest revision for a specific project.

```console
$ nextflow pull nextflow-io/hello -r v1.1

Checking nextflow-io/hello ...
checkout-out at AnyObjectId[1c3e9e7404127514d69369cd87f8036830f5cf64] - revision: 1c3e9e7404 [v1.1]
```

(cli-run)=

### run

Execute a pipeline.

**Usage**

```console
$ nextflow run [options] [project]
```

**Description**

The `run` command is used to execute a local pipeline script or remote pipeline project.

**Options**

`-E`
: Exports all current system environment.

`-ansi-log`
: Enable/disable ANSI console logging.

`-bucket-dir`
: Remote bucket where intermediate result files are stored. When running a hybrid workflow, `-bucket-dir` and `-work-dir` should define separate work directories for remote tasks and local tasks, respectively.

`-cache`
: Enable/disable processes caching.

`-d, -deep`
: Create a shallow clone of the specified depth.

`-disable-jobs-cancellation`
: Prevent the cancellation of child jobs on execution termination

`-dsl1`
: :::{deprecated} 23.09.0-edge
  :::
: Execute the workflow using DSL1 syntax.

`-dsl2`
: :::{deprecated} 23.09.0-edge
  :::
: Execute the workflow using DSL2 syntax.

`-dump-channels`
: Dump channels for debugging purpose.

`-dump-hashes`
: Dump task hash keys for debugging purposes.
: :::{versionadded} 23.10.0
  You can use `-dump-hashes json` to dump the task hash keys as JSON for easier post-processing. See the {ref}`caching and resuming tips <cache-compare-hashes>` for more details.
  :::

`-e.<key>=<value>`
: Add the specified variable to execution environment.

`-entry`
: Entry workflow to be executed.

`-h, -help`
: Print the command usage.

`-hub` (`github`)
: Service hub where the project is hosted. Options: `gitlab` or `bitbucket`

`-latest`
: Pull latest changes before run.

`-lib`
: Library extension path.

`-main-script` (`main.nf`)
: :::{versionadded} 20.09.1-edge
  :::
: The script file to be executed when launching a project directory or repository.

`-name`
: Assign a mnemonic name to the a pipeline run.

`-offline`
: Do not check for remote project updates.

`-params-file`
: Load script parameters from a JSON/YAML file.

`-plugins`
: Comma separated list of plugin ids to be applied in the pipeline execution.

`-preview`
: :::{versionadded} 22.06.0-edge
  :::
: Run the workflow script skipping the execution of all processes.

`-process.<key>=<value>`
: Set process config options.

`-profile`
: Choose a configuration profile.

`-qs, -queue-size`
: Max number of processes that can be executed in parallel by each executor.

`-resume`
: Execute the script using the cached results, useful to continue executions that was stopped by an error.

`-r, -revision`
: Revision of the project to run (either a git branch, tag or commit SHA number).
: When passing a git tag or branch, the `workflow.revision` and `workflow.commitId` fields are populated. When passing only the commit hash, `workflow.revision` is not defined.

`-stub-run, -stub`
: Execute the workflow replacing process scripts with command stubs

`-test`
: Test a script function with the name specified.

`-user`
: Private repository user name.

`-with-apptainer`
: Enable process execution in an Apptainer container.

`-with-charliecloud`
: Enable process execution in a Charliecloud container.

`-with-cloudcache`
: Enable the use of the Cloud cache plugin for storing cache metadata to an object storage bucket.

`-with-conda`
: Use the specified Conda environment package or file (must end with `.yml` or `.yaml`)

`-with-dag` (`dag-<timestamp>.html`)
: Create pipeline DAG file.
: :::{versionchanged} 23.10.0
  The default format was changed from `dot` to `html`.
  :::

`-with-docker`
: Enable process execution in a Docker container.

`-N, -with-notification`
: Send a notification email on workflow completion to the specified recipients.

`-with-podman`
: Enable process execution in a Podman container.

`-with-report` (`report-<timestamp>.html`)
: Create workflow execution HTML report.

`-with-singularity`
: Enable process execution in a Singularity container.

`-with-spack`
: Use the specified Spack environment package or file (must end with `.yaml`)

`-with-timeline` (`timeline-<timestamp>.html`)
: Create workflow execution timeline.

`-with-tower` (`https://api.cloud.seqera.io`)
: Monitor workflow execution with [Seqera Platform](https://seqera.io/) (formerly Tower Cloud).

`-with-trace` (`trace-<timestamp>.txt`)
: Create workflow execution trace file.

`-with-wave` (`https://wave.seqera.io`)
: Enable the use of Wave containers.

`-with-weblog` (`http://localhost`)
: Send workflow status messages via HTTP to target URL.

`-without-conda`
: Disable process execution with Conda.

`-without-docker`
: Disable process execution with Docker.

`-without-podman`
: Disable process execution in a Podman container.

`-without-spack`
: Disable process execution with Spack.

`-without-wave`
: Disable the use of Wave containers.

`-w, -work-dir` (`work`)
: Directory where intermediate result files are stored.

**Examples**

- Run a specific revision of a remote pipeline.

  ```console
  $ nextflow run nextflow-io/hello -r v1.1

  N E X T F L O W  ~  version 20.07.1
  Launching `nextflow-io/hello` [grave_cajal] - revision: 1c3e9e7404 [v1.1]
  ```

- Choose a `profile` for running the project. Assumes that a profile named `docker` has already been defined in the config file.

  ```console
  $ nextflow run main.nf -profile docker
  ```

- Execute a pipeline and generate the summary HTML report. For more information on the metrics, please refer the {ref}`tracing-page` section:

  ```console
  $ nextflow run main.nf -with-report
  ```

- Execute a pipeline with a custom queue size. By default, the queue size is the number of available CPUs.

  ```console
  $ nextflow run nextflow-io/hello -qs 4
  ```

- Invoke the pipeline with a specific workflow as the entry-point.

  ```console
  $ nextflow run main.nf -entry workflow_A
  ```

- Execute a pipeline with integrated monitoring in [Seqera Platform](https://seqera.io).

  ```console
  $ nextflow run nextflow-io/hello -with-tower
  ```

- Execute a pipeline with a custom parameters file (YAML or JSON).

  ```console
  $ nextflow run main.nf -params-file pipeline_params.yml
  ```

  For example, the following params file in YAML format:

  ```yaml
  alpha: 1
  beta: 'foo'
  ```

  Or in JSON format:

  ```json
  {
    "alpha": 1,
    "beta": "foo"
  }
  ```

  Is equivalent to the following command line:

  ```console
  $ nextflow run main.nf --alpha 1 --beta foo
  ```

  The parameters specified with this mechanism are merged with the resolved configuration (base configuration and profiles). The values provided via a params file overwrite those of the same name in the Nextflow configuration file.

### self-update

Update the nextflow runtime to the latest available version.

**Usage**

```console
$ nextflow self-update
```

**Description**

The `self-update` command directs the `nextflow` CLI to update itself to the latest stable release.

**Examples**

Update Nextflow.

```console
$ nextflow self-update

      N E X T F L O W
      version 20.07.1 build 5412
      created 24-07-2020 15:18 UTC (20:48 IDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io


Nextflow installation completed. Please note:
- the executable file `nextflow` has been created in the folder: /usr/local/bin
```

### view

View a project's script file(s).

**Usage**

```console
$ nextflow view [options] [project]
```

**Description**

The `view` command is used to inspect the pipelines that are already stored in the global nextflow cache. For downloading a pipeline into the global cache `~/.nextflow/assets`, refer to the `pull` command.

**Options**

`-h, -help`
: Print the command usage.

`-l`
: List repository content.

`-q`
: Hide header line.

`-r, -revision`
: Revision of the project (either a git branch, tag or commit SHA number).

**Examples**

Viewing the contents of a downloaded pipeline.

```console
$ nextflow view nextflow-io/hello

== content of file: .nextflow/assets/nextflow-io/hello/main.nf
#!/usr/bin/env nextflow

process sayHello {
  input:
    val x
  output:
    stdout
  script:
    """
    echo '$x world!'
    """
}

workflow {
  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola') | sayHello | view
}
```

List the folder structure of the downloaded pipeline:

```console
$ nextflow view -l nextflow-io/hello

== content of path: .nextflow/assets/nextflow-io/hello
LICENSE
README.md
nextflow.config
.gitignore
circle.yml
foo.nf
.git
.travis.yml
main.nf
```

View the contents of a downloaded pipeline without omitting the header:

```console
$ nextflow view -q nextflow-io/hello

#!/usr/bin/env nextflow

process sayHello {
  input:
    val x
  output:
    stdout
  script:
    """
    echo '$x world!'
    """
}

workflow {
  Channel.of('Bonjour', 'Ciao', 'Hello', 'Hola') | sayHello | view
}
```

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
