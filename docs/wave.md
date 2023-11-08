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
The Tower access token is not mandatory, but it is recommended in order to access private container repositories and pull public containers without being affected by service rate limits. Credentials should be made available to Wave using the [credentials manager](https://help.tower.nf/latest/credentials/overview) in Tower.
:::

## Use cases

(wave-authenticate-private-repos)=

### Authenticate private repositories

Wave allows the use of private repositories in your Nextflow pipelines. The repository access keys must be provided in the form of [Tower credentials](https://help.tower.nf/latest/credentials/overview/).

Once the credentials have been created, simply specify your [Tower account access token](https://help.tower.nf/22.2/api/overview/#authentication) in your pipeline configuration file. If the credentials were created in a Tower organization workspace, specify the workspace ID as well in the config file as shown below:

```groovy
tower {
  accessToken = '<your access token>'
  workspaceId = '<your workspace id>'
}
```

### Build module containers

Wave can build and provision container images on-demand for your Nextflow pipelines.

To enable this feature, add the Dockerfile of the container to be built in the {ref}`module directory <module-directory>` where the pipeline process is defined. When Wave is enabled, it automatically uses the Dockerfile to build the required container, upload to the registry, and it uses the container to carry out the tasks defined in the module.

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

### Build Singularity native images

:::{versionadded} 23.09.0-edge
:::

Nextflow can build Singularity native images on-demand either using `Singularityfile`,
Conda packages or Spack packages. The Singularity images are automatically uploaded in a container registry OCI compliant
of your choice and stored as a [ORAS artefact](https://oras.land/).

:::{note}
This feature requires of Singularity (or Apptainer) version supporting the pull of images using the `oras:` pseudo-protocol.
:::

For example to enable the provisioning of Singularity images in your pipeline use the following configuration snippet:

```groovy
singularity.enabled = true
wave.enabled = true
wave.freeze = true
wave.strategy = ['conda']
wave.build.repository = 'docker.io/user/repo'
```

In the above configuration replace `docker.io/user/repo` with a repository of your choice where Singularity image files
should be uploaded.

:::{note}
When using a private repository, the repository access keys must be provider via Tower credentials manager (see {ref}`above <wave-authenticate-private-repos>`).

Moreover the access to the repository must be granted in the compute nodes by using the command `singularity remote login <registry>`.
Please see Singularity documentation for further details.
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

`wave.freeze`
: :::{versionadded} 23.07.0-edge
  :::
: When enabling the container freeze mode, Wave will provision an non-ephemeral container image
that will be pushed to a container repository your choice. It requires the use of the `wave.build.repository` setting.
It is also suggested to specify a custom cache repository via the setting `wave.build.cacheRepository`. Note: when using
container freeze mode, the container repository authentication needs to be managed by the underlying infrastructure.    

`wave.build.repository`
: The container repository where images built by Wave are uploaded (note: the corresponding credentials must be provided in your Nextflow Tower account).

`wave.build.cacheRepository`
: The container repository used to cache image layers built by the Wave service (note: the corresponding credentials must be provided in your Nextflow Tower account).

`wave.build.conda.basePackages`
: One or more Conda packages to be always added in the resulting container (default: `conda-forge::procps-ng`).

`wave.build.conda.commands`
: One or more commands to be added to the Dockerfile used to build a Conda based image.

`wave.build.conda.mambaImage`
: The Mamba container image is used to build Conda based container. This is expected to be [micromamba-docker](https://github.com/mamba-org/micromamba-docker) image.

`wave.build.spack.basePackages`
: :::{versionadded} 22.06.0-edge
  :::
: One or more Spack packages to be always added in the resulting container.

`wave.build.spack.commands`
: :::{versionadded} 22.06.0-edge
  :::
: One or more commands to be added to the Dockerfile used to build a Spack based image.

`wave.httpClient.connectTime`
: :::{versionadded} 22.06.0-edge
  :::
: Sets the connection timeout duration for the HTTP client connecting to the Wave service (default: `30s`).

`wave.strategy`
: The strategy to be used when resolving ambiguous Wave container requirements (default: `'container,dockerfile,conda,spack'`).

`wave.report.enabled` (preview)
: :::{versionadded} 23.06.0-edge
  :::
: Enable the reporting of the Wave containers used during the pipeline execution (default: `false`).

`wave.report.file` (preview)
: :::{versionadded} 23.06.0-edge
  :::
: The name of the containers report file (default: `'containers-<timestamp>.config'`).

`wave.retryPolicy.delay`
: :::{versionadded} 22.06.0-edge
  :::
: The initial delay when a failing HTTP request is retried (default: `150ms`). 

`wave.retryPolicy.maxDelay`
: :::{versionadded} 22.06.0-edge
  :::
: The max delay when a failing HTTP request is retried (default: `90 seconds`).

`wave.retryPolicy.maxAttempts`
: :::{versionadded} 22.06.0-edge
  :::
: The max number of attempts a failing HTTP request is retried (default: `5`).

`wave.retryPolicy.jitter`
: :::{versionadded} 22.06.0-edge
  :::
: Sets the jitterFactor to randomly vary retry delays by (default: `0.25`).
