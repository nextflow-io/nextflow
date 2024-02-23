(sharing-page)=

# Sharing pipelines

Nextflow seamlessly integrates with popular Git providers, including [BitBucket](http://bitbucket.org/), [GitHub](http://github.com), and [GitLab](http://gitlab.com) for sharing pipelines. This feature allows you to version control your pipeline code, use other people's Nextflow pipelines, and publish your pipelines, in a quick and transparent way.

## How it works

When you launch a script execution with Nextflow, it will look for a file with the pipeline name you've specified. If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified). If it is found, the repository is automatically downloaded to your computer and executed. This repository is stored in the Nextflow home directory, that is by default the `$HOME/.nextflow` path, and thus will be reused for any further executions.

## Running a pipeline

To launch the execution of a pipeline project, hosted in a remote code repository, you simply need to specify its qualified name or the repository URL after the `run` command. The qualified name is formed by two parts: the `owner` name and the `repository` name separated by a `/` character.

In other words if a Nextflow project is hosted, for example, in a GitHub repository at the address `http://github.com/foo/bar`, it can be executed by entering the following command in your shell terminal:

```bash
nextflow run foo/bar
```

or using the project URL:

```bash
nextflow run http://github.com/foo/bar
```

:::{note}
In the first case, if your project is hosted on a service other than GitHub, you will need to specify this hosting service in the command line by using the `-hub` option. For example `-hub bitbucket` or `-hub gitlab`. In the second case, i.e. when using the project URL as name, the `-hub` option is not needed.
:::

You can try this feature out by simply entering the following command in your shell terminal:

```bash
nextflow run nextflow-io/hello
```

It will download a trivial `Hello` example from the repository published at the following address <http://github.com/nextflow-io/hello> and execute it in your computer.

If the `owner` part in the pipeline name is omitted, Nextflow will look for a pipeline between the ones you have already executed having a name that matches the name specified. If none is found it will try to download it using the `organisation` name defined by the environment variable `NXF_ORG` (which by default is `nextflow-io`).

:::{tip}
To access a private repository, specify the access credentials by using the `-user` command line option, then the program will ask you to enter the password interactively. Private repository access credentials can also be defined in the [Git configuration file](#git-configuration).
:::

## Handling revisions

Any Git branch, tag or commit ID defined in a project repository, can be used to specify the revision that you want to execute when launching a pipeline by adding the `-r` option to the run command line. So for example you could enter:

```bash
nextflow run nextflow-io/hello -r mybranch
```

or

```bash
nextflow run nextflow-io/hello -r v1.1
```

It will execute two different project revisions corresponding to the Git tag/branch having that names.

## Commands to manage projects

The following commands allows you to perform some basic operations that can be used to manage your projects.

:::{note}
Nextflow is not meant to completely replace the [Git](https://git-scm.com/) tool. You may still need `git` to create new repositories or commit changes, etc.
:::

### Listing available projects

The `list` command allows you to list all the projects you have downloaded in your computer. For example:

```bash
nextflow list
```

This prints a list similar to the following one:

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

Starting from the top it shows: 1) the project name; 2) the Git repository URL; 3) the local directory where the project has been downloaded; 4) the script that is executed when launched; 5) the list of available revisions i.e. branches and tags. Tags are marked with a `[t]` on the right, the current checked-out revision is marked with a `*` on the left.

### Pulling or updating a project

The `pull` command allows you to download a project from a GitHub repository or to update it if that repository has already been downloaded. For example:

```bash
nextflow pull nextflow-io/examples
```

Alternatively, you can use the repository URL as the name of the project to pull:

```bash
nextflow pull https://github.com/nextflow-io/examples
```

Downloaded pipeline projects are stored in the directory `$HOME/.nextflow/assets` in your computer.

### Viewing the project code

The `view` command allows you to quickly show the content of the pipeline script you have downloaded. For example:

```bash
nextflow view nextflow-io/hello
```

By adding the `-l` option to the example above it will list the content of the repository.

### Cloning a project into a directory

The `clone` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example:

```bash
nextflow clone nextflow-io/hello target-dir
```

If the destination directory is omitted the specified project is cloned to a directory with the same name as the pipeline base name (e.g. `hello`) in the current directory.

The clone command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push back your changes by using the usual Git/GitHub workflow.

### Deleting a downloaded project

Downloaded pipelines can be deleted by using the `drop` command, as shown below:

```bash
nextflow drop nextflow-io/hello
```

(sharing-scm-file)=

## Git configuration

You can configure your credentials for various Git providers in the SCM configuration file, located at `$HOME/.nextflow/scm`. Refer to the {ref}`git-page` page for more information.

## Local repository configuration

Nextflow is also able to handle repositories stored in a local or shared file system. The repository must be created as a [bare repository](https://mijingo.com/blog/what-is-a-bare-git-repository).

Having, for example. a bare repository store at path `/shared/projects/foo.git`, Nextflow is able to run it using the following syntax:

```bash
nextflow run file:/shared/projects/foo.git
```

See [Git documentation](https://git-scm.com/book/en/v2/Git-on-the-Server-Getting-Git-on-a-Server) for more details about how create and manage bare repositories.

## Publishing your pipeline

In order to publish your Nextflow pipeline to GitHub (or any other supported platform) and allow other people to use it, you only need to create a GitHub repository containing all your project script and data files. If you don't know how to do it, follow this simple tutorial that explains how [create a GitHub repository](https://help.github.com/articles/create-a-repo).

Nextflow only requires that the main script in your pipeline project is called `main.nf`. A different name can be used by specifying the `manifest.mainScript` attribute in the `nextflow.config` file that must be included in your project. For example:

```groovy
manifest.mainScript = 'my_very_long_script_name.nf'
```

To learn more about this and other project metadata information, that can be defined in the Nextflow configuration file, read the {ref}`Manifest <config-manifest>` section on the Nextflow configuration page.

Once you have uploaded your pipeline project to GitHub other people can execute it simply using the project name or the repository URL.

For if your GitHub account name is `foo` and you have uploaded a project into a repository named `bar` the repository URL will be `http://github.com/foo/bar` and people will able to download and run it by using either the command:

```bash
nextflow run foo/bar
```

or

```bash
nextflow run http://github.com/foo/bar
```

See the [Running a pipeline](#running-a-pipeline) section for more details on how to run Nextflow projects.

## Managing dependencies

Computational pipelines are rarely composed by a single script. In real world applications they depend on many other components, including other scripts and tools, databases, and specialized environments which provide compute and storage.

These external dependencies are the primary challenge when sharing software, because the users need to recreate the environment around a tool in order to use it. This setup process is often painful and error prone, which severely hinders the ability to reproduce computational results on a system other than the one on which it was originally developed.

Nextflow tackles this problem by integrating with existing tools for reproducible software, namely Git for source code and Docker for containers. These tools allow you to keep all the dependencies of your pipeline project in one place and track changes over time with version control.

By making your pipeline project is self-contained, meaning all of its dependencies are fully defined in the project itself, you gain two major advantages:

- **Portability**: the pipeline can be run in virtually any environment with a Java VM and a container runtime
- **Reproducibility**: any results produced by the pipelined can be easily reproduced, even across different environments

One way to account for dependencies is to break them down into three categories: code, data, and environment. Here we will describe how to include each of these dependencies in your Nextflow pipeline:

### Code

Aside from pipeline scripts, you may have additional scripts and tools used by individual tasks.

#### Standard dependencies

Many standard tools can be accessed as Docker containers or Conda/Spack packages. In this case, you only need to specify the container image URL (e.g. from [DockerHub](https://hub.docker.com)) or package name with the `container` or `conda` directive, respectively.

Make sure to enable the desired method in your Nextflow configuration:

```groovy
// containers
docker.enabled = true

// conda packages
conda.enabled = true
```

This way, when you launch your pipeline, Nextflow will automatically download the necessary dependencies to run your tasks based on this configuration.

Read the {ref}`container-page` page to learn more about how to use containers with Nextflow, and the {ref}`conda-page` page for Conda packages.

:::{tip}
For maximal reproducibility, make sure to define a specific version for each tool. Otherwise, your pipeline might use different versions across subsequent runs, which can introduce subtle differences to your results.
:::

(bundling-executables)=

#### The `bin` directory

As for custom scripts, you can include executable scripts in the `bin` directory of your pipeline repository. When configured correctly, these scripts can be executed like a regular command from any process script (i.e. without modifying the `PATH` environment variable or using an absolute path), and changing the script will cause the task to be re-executed on a resumed run (i.e. just like changing the process script itself).

To configure a custom script:

1. Save the script in the `bin` directory (relative to the pipeline repository root).
2. Specify a portable shebang (see note below for details).
3. Make the script executable. For example: `chmod a+x bin/my_script.py`

:::{tip}
To maximize the portability of your bundled script, use `env` to dynamically resolve the location of the interpreter instead of hard-coding it in the shebang line.

For example, shebang definitions `#!/usr/bin/python` and `#!/usr/local/bin/python` both hard-code specific paths to the Python interpreter. Instead, the following approach is more portable:

```bash
#!/usr/bin/env python
```
:::

#### The `lib` directory

Any Groovy scripts or JAR files in the `lib` directory will be automatically loaded and made available to your pipeline scripts. The `lib` directory is a useful way to provide utility code or external libraries without cluttering the pipeline scripts.

### Data

In general, input data should be provided by external sources using parameters which can be controlled by the user. This way, a pipeline can be easily reused to process different datasets which are appropriate for the pipeline.

Parameters can be declared with default values in the main script or in the configuration file:

```groovy
params.my_input = 'default input file'
params.my_output = 'default output path'
params.my_flag = false
// ...
```

When launching a pipeline, parameter values can be provided on the command line or in a params file (using the `-params-file` option). Options prefixed with a double dash (`--`) are interpreted as parameters:

```bash
nextflow run <your pipeline> --my_input /path/to/input/file --my_output /other/path --my_flag true
```

In some cases, when a pipeline requires some small data which rarely changes, it may be easier to include this data in the pipeline repository. You can reference this data from the pipeline script in a portable manner (i.e. without relying on an absolute path) by using the `projectDir` implicit variable, which refers to the local copy of the pipeline repository.

The following example references the file `dataset/sequences.fa` in the pipeline repository:

```groovy
sequences = file("$projectDir/dataset/sequences.fa")
sequences.splitFasta {
    println it
}
```

### Environment

The "environment" refers to any other aspects of the environment in which your pipeline is executed, such as environment variables and resource managers like SLURM.

Any environment variable that may be required by the tools in your pipeline can be defined under the `env` scope in your Nextflow configuration. For example:

```groovy
env {
  DELTA = 'foo'
  GAMMA = 'bar'
}
```

Similarly, if you use an HPC scheduler like SLURM or a cloud batch service like AWS Batch to execute tasks in a distributed manner, you can use a configuration profile to define the settings for a given environment.

See the {ref}`config-page` page to learn more about Nextflow configuration, and the {ref}`executor-page` page to learn more about specific executors.
