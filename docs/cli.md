(cli-page)=

# Command line interface

Nextflow provides a robust command line interface (CLI) for the management and execution of pipelines. This page explains the key concepts and common usage patterns for the CLI.

For the complete reference of all commands, subcommands, and options, see {ref}`cli-reference`.

:::{note}
Nextflow uses two types of command line flags:
- Nextflow options use a single dash (e.g., `-log`) and modify Nextflow's behavior.
- Pipeline parameters use a double dash (e.g., `--input`) and are passed to your pipeline script.
:::

## Pipeline execution

Pipeline execution is the core function of Nextflow. These commands run Nextflow workflows, either from local files or remote Git repositories. Nextflow handles downloading, caching, and executing pipelines with minimal user intervention.

### Launching a project

The `run` command executes pipeline scripts from local files or remote repositories. It automatically manages repository downloads, caching, and execution, supporting various Git providers and authentication methods.

See {ref}`cli-run` for more information.

**Local pipelines**

Run a pipeline from your local filesystem:

```console
$ nextflow run main.nf
```

**Remote pipelines**

Use the format `<organization>/<repository>` to run a pipeline directly from Git repositories:

```console
$ nextflow run nextflow-io/hello
```

Nextflow automatically:

1. Downloads the repository to `$HOME/.nextflow/assets/`
2. Caches it for future runs
3. Executes the main script

If you omit the organization, Nextflow searches cached pipelines first, then attempts to download from the `NXF_ORG` organization (default: `nextflow-io`).

You can also use full repository URLs:

```console
$ nextflow run https://github.com/nextflow-io/hello
```

**Private repositories**

Use the `-user` option to add credentials for private repositories:

```console
$ nextflow run organization/private-repo -user my-username
```

Alternatively, configure Git authentication. See {ref}`Git configuration <git-page>` for more information.

**Non-GitHub providers**

Use the `-hub` option specify Bitbucket, GitLab, or other Git providers:

```console
$ nextflow run organization/repo -hub bitbucket
```

**Revision selection**

Use the `-r` option to specify Git branches, tags, or commits:

```console
$ nextflow run nextflow-io/hello -r v1.1
$ nextflow run nextflow-io/hello -r dev-branch
$ nextflow run nextflow-io/hello -r a3f5c8e
```

(cli-params)=

### Pipeline parameters

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

For complex parameter sets, use YAML or JSON files with `-params-file`. This is cleaner than long command lines.

```json
{
  "input": "data.csv",
  "output": "results/",
  "min_quality": 20
}
```

```console
$ nextflow run main.nf -params-file params.json
```

**Parameter precedence**

Nextflow applies parameters defined in multiple places in the following order (lowest to highest priority):

1. Script parameters (`params.foo = 'default'`)
2. Configuration parameters (see {ref}`config-params`)
3. Parameter files (`-params-file`)
4. Command line parameters (`--foo bar`)

### Kubernetes execution

The `kuberun` command executes pipelines entirely within a Kubernetes cluster. This experimental feature runs the Nextflow driver itself inside Kubernetes.

Use this when you want both the Nextflow driver and tasks running in Kubernetes. This differs from using the Kubernetes executor with `run`, where only tasks run in Kubernetes while the driver runs externally.

```console
$ nextflow kuberun nextflow-io/hello
```

See {ref}`cli-kuberun` for more information.

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

Use this to manually download pipelines before running them, update cached pipelines, or download pipelines for offline use.

```console
$ nextflow pull nextflow-io/hello
```

You can specify a particular revision to download.

```console
$ nextflow pull nextflow-io/hello -r mybranch
```

See {ref}`cli-pull` for more information.

### Viewing project code

The `view` command displays the contents of a pipeline's main script or lists all files in the repository.

Use this to quickly inspect pipeline code without opening files or explore the project structure. Specify `-l` option lists all repository files instead of showing script contents.

```console
$ nextflow view nextflow-io/hello
$ nextflow view nextflow-io/hello -l
```

See {ref}`cli-view` for more information.

### Cloning projects locally

The `clone` command copies a pipeline from the cache to a local directory and creates a full Git repository you can modify.

Use this when you want to modify an existing pipeline, create a derivative pipeline, or study a pipeline's structure. If you omit the target directory it uses the pipeline name.

```console
$ nextflow clone nextflow-io/hello my-hello
```

See {ref}`cli-clone` for more information.

### Deleting cached projects

The `drop` command removes a downloaded pipeline from the local cache.

Use this to free disk space by removing pipelines you no longer need. The project is deleted from `$HOME/.nextflow/assets/`. The next `run` will download it again if needed.

```console
$ nextflow drop nextflow-io/hello
```

See {ref}`cli-drop` for more information.

### Secret management

The `secrets` command manages secure pipeline secrets.

Use this to store credentials securely, reference them in pipelines without exposing values, and manage sensitive data centrally across your organization.

```console
$ nextflow secrets list
$ nextflow secrets set AWS_ACCESS_KEY_ID
$ nextflow secrets delete AWS_ACCESS_KEY_ID
```

See {ref}`cli-secrets` for more information. 

## Configuration and validation

Configuration and validation options and commands help you control and verify pipeline settings. Configuration options supplement pipeline configuration at runtime, while validation commands inspect how Nextflow interprets your configuration files, process definitions, and scripts.

Use these to customize pipeline configuration, debug configuration issues, verify settings, and catch issues before execution.

### Soft configuration override

The `-c` option adds your configuration on top of the defaults and merges them together.

Use this when you want to override specific settings while keeping other defaults intact. Multiple configuration files can be specified as a comma separated list.

```console
$ nextflow -c my.config run nextflow-io/hello
```

See {ref}`config-page` for more information.

### Hard configuration override

The `-C` option replaces all default configuration with your custom configuration files. Multiple configuration files can be specified as a comma separated list.

Use this when you want to ensure no default configurations interfere with your custom settings. Unlike `-c` which merges configurations, `-C` ensures only your specified file is used.

```console
$ nextflow -C my.config run nextflow-io/hello
```

See {ref}`config-page` for more information.

### Configuration inspection

The `config` command prints the resolved configuration for a pipeline.

Use this to debug configuration issues, verify which settings will be applied, understand configuration precedence, or inspect specific configuration properties.

```console
$ nextflow config
$ nextflow config nextflow-io/hello
```

See {ref}`cli-config` for more information.

### Process inspection

:::{versionadded} 23.10.0
:::

The `inspect` command analyzes process settings in a pipeline without executing it. It outputs container information in JSON or Nextflow configuration format.

Use this to determine which container images will be used by each process before running the pipeline.

```console
$ nextflow inspect nextflow-io/hello
$ nextflow inspect nextflow-io/hello -format json
```

See {ref}`cli-inspect` for more information.

### Script validation

The `lint` command analyzes Nextflow scripts and configuration files for syntax errors and code issues. It can also automatically format your code to maintain consistent style across your project.

Use this to catch syntax errors before execution, enforce consistent code formatting, or validate entire directories of Nextflow code.

```console
$ nextflow lint main.nf
```

See {ref}`cli-lint` for more information.

## Execution history

Execution history and maintenance commands manage past runs and clean up cached files. Nextflow maintains metadata about all executions and stores intermediate files in work directories.

Use these commands to review past executions, free disk space, troubleshoot failures, or explore data lineage.

### Execution logs

The `log` command displays execution history and details about past pipeline runs, such as run names, timestamps, and customizable output fields.

Use this to find run names for resuming, review execution history, or debug failed runs. The command shows recent executions by default, with options to view specific runs or customize output fields.

```console
$ nextflow log
$ nextflow log dreamy_euler
$ nextflow log last -f name,status,duration
```

See {ref}`cli-log` for more information.

### Work directory cleanup

The `clean` command removes work directories and cached intermediate files from past executions.

Use this to free disk space, clean up failed or test runs, or maintain your work directory. Use `-n` to perform a dry run and show what would be deleted. Use `-f` to delete files.

```console
$ nextflow clean -n
$ nextflow clean dreamy_euler -f
```

See {ref}`cli-clean` for more information.

### Data lineage

:::{versionadded} 25.04.0
:::

:::{warning} *Experimental: may change in a future release.*
:::

The `lineage` command explores data lineage and provenance for workflow executions and tracks relationships between inputs, outputs, and processing steps.

Use this to understand input/output relationships between tasks, trace data flow through the pipeline, or establish file provenance. Lineage tracking must be enabled in configuration.

```console
$ nextflow lineage
```

See {ref}`data-lineage-page` to get started and {ref}`cli-lineage` for more information. 

## Seqera Platform

[Seqera Platform](https://seqera.io) is a comprehensive workflow orchestration platform that extends Nextflow with features for workflow management, monitoring, and collaboration.

Use these commands to authenticate with Seqera Platform and launch workflows directly to the Platform's managed infrastructure.

### Platform authentication

:::{versionadded} 25.10.0
:::

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

See {ref}`cli-auth` for more information.

### Platform workflow submission

:::{versionadded} 25.10.0
:::

The `launch` command submits a workflow to run on Seqera Platform's infrastructure instead of your local machine.

Use this to leverage Platform's cloud resources, monitoring capabilities, and execution management. The Platform handles resource provisioning, execution monitoring, and result storage.

```console
$ nextflow launch nextflow-io/hello
```

See {ref}`cli-launch` for more information.

## System utilities

System utilities provide administrative and development tools for managing Nextflow itself, interacting with remote filesystems, working with plugins, and debugging.

Use these commands for system administration, development, and testing.

### Interactive console

The `console` command launches an interactive Groovy console with Nextflow's execution context loaded.

Use this to test Nextflow DSL code interactively, debug expressions, explore Nextflow's APIs, or experiment with syntax. It opens a GUI or REPL depending on your environment.

```console
$ nextflow console
```

### Remote filesystem operations

The `fs` command performs filesystem operations on remote storage systems supported by Nextflow, such as S3, Google Cloud Storage, and Azure Blob Storage.

Use this to manage remote data, test cloud storage access, or perform bulk file operations without additional tools.

```console
$ nextflow fs list s3://my-bucket/data/
$ nextflow fs cat s3://my-bucket/data/file.txt
$ nextflow fs cp s3://my-bucket/data/file.txt s3://dest/
$ nextflow fs delete s3://my-bucket/data/file.txt
```

See {ref}`cli-fs` for more information.

### Plugins

The `plugin` command creates plugins, installs them, and executes plugin-specific operations.

See {ref}`cli-plugin` for more information.

**Plugin creation**

:::{versionadded} 25.04.0
:::

Use the `create` subcommand to create a new plugin scaffold for development:

```console
$ nextflow plugin create
```

**Plugin installation**

Use the `install` subcommand to install a plugin and extend Nextflow functionality:

```console
$ nextflow plugin install my-plugin
```

**Plugin execution**

Use the the format `plugin-name:command` to execute plugin-specific commands:

```console
$ nextflow plugin my-plugin:hello --alpha --beta
```

See individual plugin documentation for plugin specific commands.

### Nextflow updates

The `self-update` command updates Nextflow to a newer version.  It downloads and installs the latest release or a specific version.

Use this to upgrade Nextflow, switch versions, or install edge releases. By default, it updates to the latest stable release. Specify a particular version or use the `NXF_EDGE` environment variable for development releases.

```console
$ nextflow self-update
$ NXF_EDGE=1 nextflow self-update
```

### Command help

The `help` command displays detailed usage information for any Nextflow command.

Use this to learn about command-specific options, refresh your memory about syntax, or discover available features for a particular command.

```console
$ nextflow help run
```

### Version information

The `-v` and `-version` options print Nextflow version information.

Use `-v` for minimal output showing version and build number.

```console
$ nextflow -v
nextflow version 24.04.0.5917
```

Use `-version` for detailed output showing creation date, citation, and website.

```console
$ nextflow -version
N E X T F L O W
version 24.04.0 build 5917
created 03-05-2024 15:07 UTC
cite doi:10.1038/nbt.3820
http://nextflow.io
```
