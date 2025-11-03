(cli-page)=

# Command line interface

The Nextflow command line interface (CLI) is how you interact with Nextflow to run pipelines, manage projects, inspect configurations, and maintain your workflow environment. This page explains the key concepts and common usage patterns for the CLI.

For a complete technical reference of all options, see {ref}`cli-reference`.

:::{note}
Nextflow uses two types of command-line flags:

- Nextflow options use a single dash (e.g., `-log`) and modify Nextflow's behavior.
- Pipeline parameters use a double dash (e.g., `--input`) and are passed to your pipeline script.
:::

## Options

Nextflow options are modifiers that control how Nextflow itself operates. Use them to configure execution behavior, manage logging, override configurations, or get information about Nextflow.

### Hard configuration override

The `-C` option replaces all default configuration with your custom configuration file, ensuring complete control over settings.

Use this when you need absolute control and want to ensure no default configurations interfere with your custom settings. Unlike `-c` which merges configurations, `-C` ensures only your specified file is used.

```console
$ nextflow -C my.config run nextflow-io/hello
```

See {ref}`config-page` for more about configuration files.

### Soft configuration override

The `-c` option adds your configuration on top of the defaults and merges them together.

Use this when you want to override specific settings while keeping other defaults intact. You can specify `-c` multiple times to layer configurations, with later files taking precedence over earlier ones.

```console
$ nextflow -c my.config run nextflow-io/hello
```

See {ref}`config-page` for more about configuration files.

### Background execution

The `-bg` option runs Nextflow as a background process and frees your terminal for other work.

Use this for long-running pipelines when you need to continue using your terminal or want to disconnect without stopping execution. It works like `nohup`, allowing you to close your terminal without stopping the pipeline. Outputs are saved to `.nextflow.log` by default.

```console
$ nextflow -bg run nextflow-io/hello
```

### Quiet mode

The `-q` option suppresses the Nextflow banner and execution progress output.

Use this to keep your terminal output clean. It hides the ASCII banner and task progress but doesn't affect pipeline `println` statements or error messages.

```console
$ nextflow -q run nextflow-io/hello
```

### Custom log file path

The `-log` option directs Nextflow logs to a specific file instead of the default `.nextflow.log`.

Use this for organizing logs by project, integrating with logging systems, or storing logs in centralized locations.

```console
$ nextflow -log /var/log/nextflow.log run nextflow-io/hello
```

### Syslog integration

The `-syslog` option sends logs to a [Syslog](https://en.wikipedia.org/wiki/Syslog) server for centralized log management.

Use this to aggregate Nextflow logs with other system logs to enable centralized monitoring and analysis.

```console
$ nextflow -syslog localhost:1234 run nextflow-io/hello
```

### JVM properties

The `-D` option sets Java system properties that control the JVM running Nextflow.

Use this for configuring memory, encoding, or other Java-level behavior. You can specify multiple properties by using `-D` multiple times.

```console
$ nextflow -Dfile.encoding=UTF-8 run nextflow-io/hello
```

See {ref}`config-env-vars` for additional JVM configuration options.

### Help display

The `-h` option shows available commands and global options.

Use this to discover available commands or refresh your memory about command syntax.

```console
$ nextflow -h
```

### Version information

The `-v` and `-version` options show Nextflow version details.

Use `-v` for brief output showing version and build number.

```console
$ nextflow -v
nextflow version 24.04.0.5917
```

Use `-version` for detailed information including creation date, citation, and website.

```console
$ nextflow -version
N E X T F L O W
version 24.04.0 build 5917
created 03-05-2024 15:07 UTC
cite doi:10.1038/nbt.3820
http://nextflow.io
```

## Pipeline execution

Pipeline execution is the core function of Nextflow. These commands run workflow scripts, either from local files or remote Git repositories. Nextflow handles downloading, caching, and executing pipelines with minimal user intervention.

### Launching a project

The `run` command executes pipeline scripts. This is the primary command for running Nextflow workflows. It handles both local scripts and remote repositories seamlessly.

Use this command whenever you need to execute a pipeline. It automatically manages repository downloads, caching, and execution, supporting various Git providers and authentication methods.

**Local execution**

Run a pipeline script from your local filesystem:

```console
$ nextflow run main.nf
```

**Remote execution**

Nextflow executes pipelines directly from Git repositories without manual downloading. Use the format `<owner>/<repository>`:

```console
$ nextflow run nextflow-io/hello
```

Nextflow automatically:

1. Downloads the repository to `$HOME/.nextflow/assets/`
2. Caches it for future runs
3. Executes the main script

If you omit the owner, Nextflow searches cached pipelines first, then attempts to download from the `NXF_ORG` organization (default: `nextflow-io`).

You can also use full repository URLs:

```console
$ nextflow run https://github.com/nextflow-io/hello
```

**Private repositories**

For private repositories, use the `-user` option to provide credentials:

```console
$ nextflow run owner/private-repo -user myusername
```

Alternatively, configure Git authentication. See {ref}`Git configuration <git-page>`.

**Non-GitHub providers**

Use the `-hub` option for Bitbucket, GitLab, or other Git providers:

```console
$ nextflow run owner/repo -hub bitbucket
```

**Revision selection**

Execute specific Git branches, tags, or commits with the `-r` option. This is essential for running specific pipeline versions, testing development branches, and ensuring reproducibility by pinning to specific commits.

```console
$ nextflow run nextflow-io/hello -r v1.1
$ nextflow run nextflow-io/hello -r dev-branch
$ nextflow run nextflow-io/hello -r a3f5c8e
```

(cli-params)=

**Pipeline parameters**

Pipeline parameters are values defined with `params` in your script. Override them on the command line using the `--` prefix to customize pipeline behavior without modifying code.

```console
$ nextflow run main.nf --input data.csv --output results
```

Parameter names support automatic conversion between kebab-case and camelCase:

```console
$ nextflow run main.nf --input-file data.csv  # Becomes params.inputFile
```

Parameters without values are set to `true`:

```console
$ nextflow run main.nf --verbose  # params.verbose = true
```

:::{warning}
Quote parameters containing wildcards to prevent shell expansion:

```console
$ nextflow run main.nf --files "*.fasta"
```
:::

**Parameter files**

For complex parameter sets, use a YAML or JSON file with `-params-file`. This is cleaner than long command lines.

```yaml
# params.yml
input: data.csv
output: results/
min_quality: 20
```

```console
$ nextflow run main.nf -params-file params.yml
```

**Parameter precedence**

 Nextflow applies parameters defined in multiple places in the following order (lowest to highest priority):

1. Pipeline script defaults (`params.foo = 'default'`)
2. Configuration files (see {ref}`config-params`)
3. Parameter files (`-params-file`)
4. Command line parameters (`--foo bar`)

### Kubernetes execution

The `kuberun` command executes pipelines entirely within a Kubernetes cluster. This experimental feature runs the Nextflow driver itself inside Kubernetes.

Use this when you want both the Nextflow driver and tasks running in Kubernetes. This differs from using the Kubernetes executor with `run`, where only tasks run in Kubernetes while the driver runs externally.

```console
$ nextflow kuberun nextflow-io/hello
```

See {ref}`cli-kuberun` for configuration details.

## Project management

Project management commands interact with Git-hosted pipelines. Nextflow integrates with Git providers (e.g., GitHub, GitLab, and Bitbucket) to treat pipelines as versioned projects and maintains a local cache in `$HOME/.nextflow/assets/`.

Use these commands to explore available pipelines, inspect their code, maintain your cache, and clone projects.

### Listing downloaded projects

The `list` command shows all pipelines currently downloaded to your local cache. This helps you track which projects are available offline and manage your cache directory. Pipelines are stored in `$HOME/.nextflow/assets/` by default.

```console
$ nextflow list
```

### Showing project information

The `info` command displays detailed metadata about a downloaded project, including its repository location, local path, main script, and available revisions.

Use this to understand a project's structure, see available versions, or verify which revision is currently checked out.

```console
$ nextflow info hello
project name: nextflow-io/hello
repository  : https://github.com/nextflow-io/hello
local path  : $HOME/.nextflow/assets/nextflow-io/hello
main script : main.nf
revisions   :
* master (default)
  mybranch
  v1.1 [t]
  v1.2 [t]
```

This shows:

- The full project name and repository URL
- Where it's cached locally
- Which script runs by default
- Available revisions (branches and tags marked with `[t]`)
- Which revision is currently checked out (marked with `*`)

### Pulling or updating projects

The `pull` command downloads a pipeline or updates an existing one to the latest version from its Git repository.

Use this to manually download pipelines before running them, update cached pipelines, or download pipelines for offline use. You can specify a particular revision to download.

```console
$ nextflow pull nextflow-io/hello
$ nextflow pull nextflow-io/hello -r mybranch
```

See {ref}`cli-pull` for configuration options.

### Viewing project code

The `view` command displays the contents of a pipeline's main script or lists all files in the repository.

Use this to quickly inspect pipeline code without opening files or explore the project structure. The `-l` option lists all repository files instead of showing script contents.

```console
$ nextflow view nextflow-io/hello
$ nextflow view nextflow-io/hello -l
```

See {ref}`cli-view` for configuration options.

### Cloning projects locally

The `clone` command copies a pipeline from the cache to a local directory and creates a full Git repository you can modify.

Use this when you want to modify an existing pipeline, create a derivative pipeline, or study a pipeline's structure in detail. If you omit the target directory it uses the pipeline name.

```console
$ nextflow clone nextflow-io/hello my-hello
```

See {ref}`cli-clone` for configuration options.

### Deleting cached projects

The `drop` command removes a downloaded pipeline from the local cache.

Use this to free disk space by removing pipelines you no longer need. The project is deleted from `$HOME/.nextflow/assets/`. The next `run` will download it again if needed.

```console
$ nextflow drop nextflow-io/hello
```

See {ref}`cli-drop` for configuration options.

## Configuration and validation

Configuration and validation commands help you understand and verify pipeline settings before execution. These commands inspect how Nextflow interprets your configuration files, process definitions, and scripts.

Use these commands during development to debug configuration issues, verify settings, and catch potential problems early.

### Configuration inspection

The `config` command prints the resolved configuration for a pipeline to show how all configuration files merge together into the final effective settings.

Use this to debug configuration issues, verify which settings will be applied, understand configuration precedence, or inspect specific configuration properties.

```console
$ nextflow config
$ nextflow config nextflow-io/hello
```

See {ref}`cli-config` for configuration options.

### Process inspection

The `inspect` command displays detailed settings for processes in a pipeline, including directives, resources, and container images.

Use this to verify process configurations before running, understand resource allocations, or check which containers will be used. This shows CPU, memory, containers, process directives, and executor settings for each process.

```console
$ nextflow inspect main.nf
```

See {ref}`cli-inspect` for configuration options.

### Script validation

The `lint` command checks Nextflow scripts and configuration files for common issues and syntax violations.

Use this during development to maintain code quality, catch syntax errors, identify deprecated features, spot style inconsistencies, and find potential bugs before execution.

```console
$ nextflow lint main.nf
```

See {ref}`cli-lint` for configuration options.

## Execution history

Execution history and maintenance commands manage past runs and clean up cached files. Nextflow maintains metadata about all executions and stores intermediate files in work directories.

Use these commands to review past executions, free disk space, troubleshoot failures, or explore data lineage.

### Execution logs

The `log` command displays execution history and details about past pipeline runs, including run names, timestamps, and customizable output fields.

Use this to find run names for resuming, review execution history, or debug failed runs. The command shows recent executions by default, with options to view specific runs or customize output fields.

```console
$ nextflow log
$ nextflow log dreamy_euler
$ nextflow log last -f name,status,duration
```

See {ref}`cli-log` for configuration options.

### Work directory cleanup

The `clean` command removes work directories and cached intermediate files from past executions.

Use this to free disk space, clean up failed or test runs, or maintain your work directory. Use `-n` to perform a dry run and show what would be deleted. Use `-f` to delete files.

```console
$ nextflow clean -n
$ nextflow clean dreamy_euler -f
```

See {ref}`cli-clean` for configuration options.

### Data lineage

The `lineage` command explores data lineage and provenance for workflow executions and tracks relationships between inputs, outputs, and processing steps.

Use this feature to understand input/output relationships between tasks, trace data flow through the pipeline, or establish file provenance. Lineage tracking must be enabled in configuration.

```console
$ nextflow lineage
```

See {ref}`data-lineage-page` to get started and {ref}`cli-lineage` for configuration options. 

### Secret management

The `secrets` command manages secure pipeline secrets.

Use this to store credentials securely, reference them in pipelines without exposing values, and manage sensitive data centrally across your organization.

```console
$ nextflow secrets list
$ nextflow secrets set AWS_ACCESS_KEY_ID
$ nextflow secrets delete AWS_ACCESS_KEY_ID
```

See {ref}`cli-secrets` for configuration options. 

## Platform integration

Platform integration commands connect Nextflow with [Seqera Platform](https://seqera.io) for workflow management, monitoring, and collaboration. The Platform provides centralized execution management, monitoring dashboards, and secure credential storage.

Use these commands when working with Seqera Platform to authenticate, launch workflows, or manage secrets.

### Platform authentication

The `auth` command manages authentication credentials for Seqera Platform, saving access tokens for API interactions.

Use this to log in or out of the Platform, establishing or removing your authentication credentials.

```console
$ nextflow auth login
```

Additional authentication operations include checking login status, viewing configuration details, and logging out:

```console
$ nextflow auth status
$ nextflow auth config
$ nextflow auth logout
```

See {ref}`cli-auth` for configuration options.

### Platform workflow submission

The `launch` command submits a workflow to run on Seqera Platform's infrastructure instead of your local machine.

Use this to leverage Platform's cloud resources, monitoring capabilities, and execution management. The Platform handles resource provisioning, execution monitoring, and result storage.

```console
$ nextflow launch nextflow-io/hello
```

See {ref}`cli-launch` for configuration options.

## System utilities

System utilities provide administrative and development tools for managing Nextflow itself, interacting with remote filesystems, working with plugins, and debugging.

Use these commands for system administration, development, and testing.

### Interactive console

The `console` command launches an interactive Groovy console with Nextflow's execution context loaded.

Use this development tool to test Nextflow DSL code interactively, debug expressions, explore Nextflow's APIs, or experiment with syntax. It opens a GUI or REPL depending on your environment.

```console
$ nextflow console
```

### Remote filesystem operations

The `fs` command performs filesystem operations on remote storage systems supported by Nextflow, including S3, Google Cloud Storage, and Azure Blob Storage.

Use this to manage remote data, test cloud storage access, or perform bulk file operations without additional tools. Operations include listing, copying, moving, and deleting files.

```console
$ nextflow fs list s3://my-bucket/data/
$ nextflow fs cp s3://source/file.txt s3://dest/
$ nextflow fs delete s3://bucket/old-data/
```

### Plugins

The `plugin` command creates plugins, installs them, and executes plugin-specific operations.

Create a new plugin scaffold for development:

```console
$ nextflow plugin create
```

Install a plugin to extend Nextflow's functionality:

```console
$ nextflow plugin install my-plugin
```

Execute plugin-specific commands using the format `plugin-name:command`:

```console
$ nextflow plugin my-plugin:hello --alpha --beta
```

See {ref}`using-plugins-page` and {ref}`cli-plugin` for more information. See individual plugin documentation for plugin specific commands.

### Nextflow updates

The `self-update` command updates Nextflow to a newer version.  It downloads and installs the latest release or a specific version.

Use this to upgrade Nextflow, switch versions, or install edge releases. By default, it updates to the latest stable release. You can specify a particular version or use the `NXF_EDGE` environment variable for development releases.

```console
$ nextflow self-update
$ NXF_EDGE=1 nextflow self-update
```

### Command help

The `help` command displays detailed usage information for any Nextflow command.

Use this to learn about command-specific options, refresh your memory about syntax, or discover available features for a particular command.

```console
$ nextflow help run
$ nextflow help clean
```