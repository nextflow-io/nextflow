(conda-page)=

# Conda environments

[Conda](https://conda.io/) is an open source package and environment management system that simplifies the installation and the configuration of complex software packages in a platform agnostic manner.

Nextflow has built-in support for Conda that allows the configuration of workflow dependencies using Conda recipes and environment files.

This allows Nextflow applications to use popular tool collections such as [Bioconda](https://bioconda.github.io) and the [Python Package index](https://pypi.org/), whilst taking advantage of the configuration flexibility provided by Nextflow.

## Prerequisites

This feature requires the Conda or [Miniconda](https://conda.io/miniconda.html) package manager to be installed on your system.

## How it works

Nextflow automatically creates and activates the Conda environment(s) given the dependencies specified by each process.

Dependencies are specified by using the {ref}`process-conda` directive, providing either the names of the required Conda packages, the path of a Conda environment yaml file, or the path of an existing Conda environment directory.

:::{note}
Conda environments are stored on the file system. By default, Nextflow instructs Conda to save the required environments in the pipeline work directory. The same environment may be created/saved multiple times across multiple executions when using different work directories.
:::

You can specify the directory where the Conda environments are stored using the `conda.cacheDir` configuration property. When using a computing cluster, make sure to use a shared file system path accessible from all compute nodes. See the {ref}`configuration page <config-conda>` for details about Conda configuration.

:::{warning}
The Conda environment feature is not supported by executors that use remote object storage as a work directory. For example, AWS Batch.
:::

### Enabling Conda environment

:::{versionadded} 22.08.0-edge
:::

The use of Conda recipes specified using the {ref}`process-conda` directive needs to be enabled explicitly in the pipeline configuration file (i.e. `nextflow.config`):

```groovy
conda.enabled = true
```

Alternatively, it can be specified by setting the variable `NXF_CONDA_ENABLED=true` in your environment or by using the `-with-conda` command line option.

### Use Conda package names

Conda package names can specified using the `conda` directive. Multiple package names can be specified by separating them with a blank space. For example:

```nextflow
process foo {
  conda 'bwa samtools multiqc'

  script:
  """
  your_command --here
  """
}
```

Using the above definition, a Conda environment that includes BWA, Samtools and MultiQC tools is created and activated when the process is executed.

The usual Conda package syntax and naming conventions can be used. The version of a package can be specified after the package name as shown here `bwa=0.7.15`.

The name of the channel where a package is located can be specified prefixing the package with the channel name as shown here `bioconda::bwa=0.7.15`.

(conda-env-files)=
### Use Conda environment files

Conda environments can also be defined using one or more Conda environment files. This is a file that lists the required packages and channels structured using the YAML format. For example:

```yaml
name: my-env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - star=2.5.4a
  - bwa=0.7.15
```

Read the Conda documentation for more details about how to create [environment files](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-file-manually).

The path of an environment file can be specified using the `conda` directive:

```nextflow
process foo {
  conda '/some/path/my-env.yaml'

  script:
  """
  your_command --here
  """
}
```

:::{warning}
The environment file name **must** have a `.yml` or `.yaml` extension or else it won't be properly recognised.
:::

(conda-pypi)=
### Python Packages from PyPI

Conda environment files can also be used to install Python packages from the [PyPI repository](https://pypi.org/), through the `pip` package manager (which must also be explicitly listed as a required package):

```yaml
name: my-env-2
channels:
  - defaults
dependencies:
  - pip
  - pip:
    - numpy
    - pandas
    - matplotlib
```

### Conda text files

It is possible to provide dependencies by listing each package name as a separate line in a plain text file. For example:

```
bioconda::star=2.5.4a
bioconda::bwa=0.7.15
bioconda::multiqc=1.4
```

:::{warning}
Like before, the extension matters. Make sure the dependencies file has a `.txt` extension.
:::

### Conda lock files

The final way to provide packages to Conda is with [Conda lock files](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#identical-conda-envs).

These are generated from existing Conda environments using the following command:

```bash
conda list --explicit > spec-file.txt
```

or if using Mamba / Micromamba:

```bash
micromamba env export --explicit > spec-file.txt
```

Conda lock files can also be downloaded from [Wave](https://seqera.io/wave/) build pages.

These files include every package and their dependencies. As such, no Conda environment resolution step is needed. This is faster and more reproducible.

The files contain package URLs and an optional md5hash for each download to confirm identity:

```
# micromamba env export --explicit
# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: linux-64
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
https://conda.anaconda.org/conda-forge/linux-64/libgomp-13.2.0-h77fa898_7.conda#abf3fec87c2563697defa759dec3d639
https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-2_gnu.tar.bz2#73aaf86a425cc6e73fcf236a5a46396d
https://conda.anaconda.org/conda-forge/linux-64/libgcc-ng-13.2.0-h77fa898_7.conda#72ec1b1b04c4d15d4204ece1ecea5978
# .. and so on
```

To use with Nextflow, simply set the `conda` directive to the lock file path.

### Use existing Conda environments

If you already have a local Conda environment, you can use it in your workflow specifying the installation directory of such environment by using the `conda` directive:

```nextflow
process foo {
  conda '/path/to/an/existing/env/directory'

  script:
  """
  your_command --here
  """
}
```

### Use Mamba to resolve packages

:::{warning} *Experimental: may change in a future release.*
:::

It is also possible to use [mamba](https://github.com/mamba-org/mamba) to speed up the creation of conda environments. For more information on how to enable this feature please refer to {ref}`Conda <config-conda>`.

## Best practices

When a `conda` directive is used in any `process` definition within the workflow script, Conda tool is required for the workflow execution.

Specifying the Conda environments in a separate configuration {ref}`profile <config-profiles>` is therefore recommended to allow the execution via a command line option and to enhance the workflow portability. For example:

```groovy
profiles {
  conda {
    process.conda = 'samtools'
  }

  docker {
    process.container = 'biocontainers/samtools'
    docker.enabled = true
  }
}
```

The above configuration snippet allows the execution either with Conda or Docker specifying `-profile conda` or `-profile docker` when running the workflow script.

## Advanced settings

Conda advanced configuration settings are described in the {ref}`Conda <config-conda>` section on the Nextflow configuration page.
