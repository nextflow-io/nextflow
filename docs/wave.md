(wave-page)=

# Wave containers

:::{versionadded} 22.10.0
:::

[Wave](https://seqera.io/wave/) is a container provisioning service integrated with Nextflow. With Wave, you can build, upload, and manage the container images required by your data analysis workflows automatically and on-demand during pipeline execution.

## Getting started

### Nextflow installation

If you have already installed Nextflow, update to the latest version using this command:

```bash
nextflow -self-update
```

If you don't have Nextflow already installed, install it with the command below:

```bash
curl get.nextflow.io | bash
```

### Wave configuration

Wave can be used in any Nextflow pipeline by adding the following snippet to your `nextflow.config` file:

```groovy
wave {
  enabled = true
}

tower {
  accessToken = '<your access token>'
}
```

:::{note}
The Tower access token is not mandatory, but it is recommended in order to access private container repositories and pull public containers without being affected by service rate limits. Credentials should be made available to Wave using the [credentials manager](https://help.tower.nf/latest/credentials/registry_credentials/) in Tower.
:::

## Use cases

### Authenticate private repositories

Wave allows the use of private repositories in your Nextflow pipelines. The repository access keys must be provided in the form of [Nextflow Tower credentials](https://help.tower.nf/22.2/credentials/overview/).

Once the credentials have been created, simply specify your [Tower account access token](https://help.tower.nf/22.2/api/overview/#authentication) in your pipeline configuration file. If the credentials were created in a Tower organization workspace, specify the workspace ID as well in the config file as shown below:

```groovy
tower {
  accessToken = '<your access token>'
  workspaceId = '<your workspace id>'
}
```

### Build module containers

Wave can build and provision container images on-demand for your Nextflow pipelines.

To enable this feature, add the Dockerfile of the container to be built in the {ref}`module directory <dsl2-module-directory>` where the pipeline process is defined. When Wave is enabled, it automatically uses the Dockerfile to build the required container, upload to the registry, and it uses the container to carry out the tasks defined in the module.

:::{tip}
Make sure the process does not declare a `container` directive, otherwise it will take precedence over the Dockerfile definition.
:::

If a process uses a `container` directive and you still want to build the container using the Dockerfile provided in the module directory, add the following setting to the pipeline config file:

```groovy
wave.strategy = ['dockerfile','container']
```

This setting instructs Wave to prioritize the module Dockerfile over process `container` directives.

:::{warning}
When building containers, Wave currently does not support `ADD`, `COPY`, or any other Dockerfile commands that access files in the host file system.
:::

### Build Conda based containers

Wave allows the provisioning of containers based on the {ref}`process-conda` directive used by the processes in your pipeline. This is a quick alternative to building Conda packages in the local computer. Moreover, this enables the use of Conda packages in your pipeline when deploying in cloud-native platforms such as AWS Batch and Kubernetes, which do not allow the (easy) use of the Conda package manager.

With Wave enabled in your pipeline, simply define the `conda` requirements in the pipeline processes, provided the same process does not also specify a `container` directive or a Dockerfile.

In the latter case, add the following setting to your pipeline configuration:

```groovy
wave.strategy = ['conda']
```

The above setting instructs Wave to use the `conda` directive to provision the pipeline containers and ignore the `container` directive and any Dockerfile(s).

### Build Spack based containers

:::{warning}
Spack based Wave containers are currently in beta testing. Functionality is still sub-optimal, due to long build times that may result in backend time-out and subsequent task failure.
:::

Wave allows the provisioning of containers based on the {ref}`process-spack` directive used by the processes in your
pipeline. This is an alternative to building Spack packages in the local computer.
Moreover, this enables to run optimised builds with almost no user intervention.

Having Wave enabled in your pipeline, there's nothing else to do other than define the `spack` requirements in
the pipeline processes provided the same process does not also specify a `container` or `conda` directive or a Dockerfile.

In the latter case, add the following setting to your pipeline configuration:

```groovy
wave.strategy = ['spack']
```

The above setting instructs Wave to only use the `spack` directive to provision the pipeline containers, ignoring the use of
the `container` directive and any Dockerfile(s).

In order to request the build of containers that are optimised for a specific CPU microarchitecture, the latter can be specified by means of the {ref}`process-arch` directive. The architecture must always be specified for processes that run on an ARM system. Otherwise, by default, Wave will build containers for the generic `x86_64` architecture family.

:::{note}
If using a Spack YAML file to provide the required packages, you should avoid editing the following sections, which are already configured by the Wave plugin: `packages`, `config`, `view` and `concretizer` (your edits may be ignored), and `compilers` (your edits will be considered, and may interfere with the setup by the Wave plugin).
:::

### Push to a private repository

Containers built by Wave are uploaded to the Wave default repository hosted on AWS ECR at `195996028523.dkr.ecr.eu-west-1.amazonaws.com/wave/build`. The images in this repository are automatically deleted 1 week after the date of their push.

If you want to store Wave containers in your own container repository use the following settings in the Nextflow configuration file:

```groovy
wave.build.repository = 'example.com/your/build-repo'
wave.build.cacheRepository = 'example.com/your/cache-repo'
```

The first repository is used to store the built container images. The second one is used to store the individual image layers for caching purposes.

The repository access keys must be provided as Tower credentials (see
[Authenticate private repositories](#authenticate-private-repositories) above).

### Run pipelines using Fusion file system

Wave containers allows you to run your containerised workflow with the {ref}`fusion-page`.

This enables the use of an object storage bucket such as AWS S3 or Google Cloud Storage as your pipeline work directory, simplifying and speeding up many operations on local, AWS Batch, Google Batch or Kubernetes executions.

See the {ref}`Fusion documentation <fusion-page>` for more details.

## Advanced settings

The following configuration options are available:

`wave.enabled`
: Enable/disable the execution of Wave containers.

`wave.endpoint`
: The Wave service endpoint (default: `https://wave.seqera.io`).

`wave.build.repository`
: The container repository where images built by Wave are uploaded (note: the corresponding credentials must be provided in your Nextflow Tower account).

`wave.build.cacheRepositor`
: The container repository used to cache image layers built by the Wave service (note: the corresponding credentials must be provided in your Nextflow Tower account).

`wave.build.conda.mambaImage`
: The Mamba container image is used to build Conda based container. This is expected to be [micromamba-docker](https://github.com/mamba-org/micromamba-docker) image.

`wave.build.conda.commands`
: One or more commands to be added to the Dockerfile used to build a Conda based image.

`wave.build.conda.basePackages`
: One or more Conda packages to be always added in the resulting container e.g. `conda-forge::procps-ng`.

`wave.build.spack.checksum`
: Enable checksum verification for source tarballs (recommended). Disable only when requesting a package version not yet encoded in the corresponding Spack recipe (default: `true`).

`wave.build.spack.builderImage`
: The Spack container image is used to build Spack based container. This is expected to be one of the [Spack-provided](https://spack.readthedocs.io/en/latest/containers.html) images.

`wave.build.spack.runnerImage`
: The OS container image is used for the production container. This is expected to match the OS of the `builderImage` above.

`wave.build.spack.osPackages`
: Additional OS packages to be installed in the production container. Note that package names may vary depending on the OS of the `runnerImage` above.

`wave.build.spack.cFlags`
: C compiler flags used during the build. Default: `-O3` for GCC compiler. Recommended: one of `-O3` (high optimisation) or `-O2` (moderate optimisation).

`wave.build.spack.cxxFlags`
: C++ compiler flags used during the build. Default: `-O3` for GCC compiler. Recommended: one of `-O3` (high optimisation) or `-O2` (moderate optimisation).

`wave.build.spack.fFlags`
: Fortran compiler flags used during the build. Default: `-O3` for GCC compiler. Recommended: one of `-O3` (high optimisation) or `-O2` (moderate optimisation).

`wave.build.spack.commands`
: One or more commands to be added to the Dockerfile used to build a Spack based image.

`wave.strategy`
: The strategy to be used when resolving ambiguous Wave container requirements (default: `'container,dockerfile,conda,spack'`).

