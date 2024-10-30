(sharing-page)=

# Sharing pipelines

Nextflow seamlessly integrates with popular Git providers, including [BitBucket](http://bitbucket.org/), [GitHub](http://github.com), and [GitLab](http://gitlab.com) for managing Nextflow pipelines as version-controlled Git repositories. This feature allows you to easily use other people's Nextflow pipelines and publish your own pipelines.

:::{note}
Nextflow is not meant to completely replace the [Git](https://git-scm.com/) tool. You may still need `git` to create new repositories or commit changes, etc.
:::

## Git configuration

You can configure your credentials for various Git providers in the Git configuration file, located at `$HOME/.nextflow/scm`. See {ref}`git-page` for more information.

## Using a local repository

Nextflow can work with repositories stored in a local or shared file system. The repository must be created as a [bare repository](https://mijingo.com/blog/what-is-a-bare-git-repository).

For example, given a bare repository at `/shared/projects/foo.git`, Nextflow is able to run it using the following syntax:

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

See the {ref}`CLI <cli-page>` page to learn how to use the Nextflow command line to run pipelines and manage pipeline projects.

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

When a pipeline requires some small data that rarely changes, it may be easier to include the data in the pipeline repository. You can reference this data from the pipeline script in a portable manner (i.e. without relying on an absolute path) by using the `projectDir` implicit variable, which refers to the local copy of the pipeline repository.

The following example references the file `dataset/sequences.fa` in the pipeline repository:

```nextflow
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

See {ref}`config-page` for more information about Nextflow configuration and {ref}`executor-page` for more information about executors.
