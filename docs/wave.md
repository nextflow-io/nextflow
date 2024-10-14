(wave-page)=

# Wave containers

:::{versionadded} 22.10.0
:::

[Wave](https://seqera.io/wave/) is a container provisioning service integrated with Nextflow. With Wave, you can build, upload, and manage the container images required by your data analysis workflows automatically and on-demand during pipeline execution.

## Getting started

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
The Seqera Platform access token is not mandatory, but it is recommended in order to access private container repositories and pull public containers without being affected by service rate limits. Credentials should be made available to Wave using the [credentials manager](https://docs.seqera.io/platform/latest/credentials/overview) in Seqera Platform.
:::

## Use cases

(wave-authenticate-private-repos)=

### Authenticate private repositories

Wave allows the use of private repositories in your Nextflow pipelines. The repository access keys must be provided in the form of [Seqera Platform credentials](https://docs.seqera.io/platform/latest/credentials/overview/).

Once the credentials have been created, simply specify your [personal access token](https://docs.seqera.io/platform/23.3.0/api/overview#authentication) in your pipeline configuration file. If the credentials were created in a Seqera Platform organization workspace, specify the workspace ID as well in the config file as shown below:

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

:::{tip}
Some configuration options in the `conda` scope are used when Wave is used to build Conda-based containers.
For example, the Conda channels and their priority can be set with `conda.channels`:

```groovy
wave.strategy = 'conda'
conda.channels = 'conda-forge,bioconda'
```
:::

(wave-singularity)=

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
When using a private repository, the repository access keys must be provided via the Seqera Platform credentials manager (see {ref}`above <wave-authenticate-private-repos>`).

Moreover the access to the repository must be granted in the compute nodes by using the command `singularity remote login <registry>`.
Please see Singularity documentation for further details.
:::

:::{note}
In order to build Singularity native images, both `singularity.ociAutoPull` and `singularity.ociMode` need to be disabled in the configuration (see the {ref}`config-singularity` section).
:::

### Push to a private repository

Containers built by Wave are uploaded to the Wave default repository hosted on AWS ECR at `195996028523.dkr.ecr.eu-west-1.amazonaws.com/wave/build`. The images in this repository are automatically deleted 1 week after the date of their push.

If you want to store Wave containers in your own container repository use the following settings in the Nextflow configuration file:

```groovy
wave.build.repository = 'example.com/your/build-repo'
wave.build.cacheRepository = 'example.com/your/cache-repo'
```

The first repository is used to store the built container images. The second one is used to store the individual image layers for caching purposes.

The repository access keys must be provided as Seqera Platform credentials (see
[Authenticate private repositories](#authenticate-private-repositories) above).

### Mirroring containers

Wave allows mirroring, i.e., copying containers used by your pipeline into a container registry of your choice. This allows the pipeline to pull containers from the target registry rather than the original registry.

Mirroring is useful to create an on-demand cache of container images that are co-located in the same region where the pipeline
is executed, and therefore optimising cost and network efficiency.

Include the following settings in your Nextflow configuration to enable this capability:

```groovy
wave.enabled = true
wave.mirror = true
wave.build.repository = '<YOUR REGISTRY>'
tower.accessToken = '<YOUR ACCESS TOKEN>'
```

In the above snippet, replace `<YOUR REGISTRY>` with a container registry of your choice. For example, `quay.io` (no prefix or suffix is needed).
The container will be copied with the same name, tag, and checksum in the specified registry. For example, if the source
container is `quay.io/biocontainers/bwa:0.7.13--1` and the build repository setting is `foo.com`, the resulting container
name is `foo.com/biocontainers/bwa:0.7.13--1`.

:::{tip}
When using a path prefix in the target registry name, it will be prepended to the resulting container name. For example,
having `quay.io/biocontainers/bwa:0.7.13--1` as source container and `foo.com/bar` as build repository, the resulting
container will be named `foo.com/bar/biocontainers/bwa:0.7.13--1`.
:::

The credentials to allow the push of  containers in the target repository need to be provided via the Seqera Platform
credentials manager. The account used for this is specified by the `tower.accessToken` in the configuration above.

### Container security scanning

Wave enables the scanning of containers used in your pipelines for security vulnerabilities.
If any issues are detected, it will trigger an execution error and provide a report.

To enable this capability add the following settings to your Nextflow configuration file:

```groovy
wave.enabled = true
wave.scan.mode = 'required'
tower.accessToken = '<YOUR ACCESS TOKEN>'
```

Nextflow will only allow the use of containers with no security
vulnerabilities when using these settings. You can define the level of accepted vulnerabilities using `wave.scan.allowedLevels`. For example:

```
wave.scan.allowedLevels = 'low,medium'
```

The above setting will allow the use of containers with *low* and *medium* vulnerabilities. Accepted values are `low`, `medium`, `high`, and `critical`. See [common vulnerabilities scoring system](https://en.wikipedia.org/wiki/Common_Vulnerability_Scoring_System) for more information about these levels.

:::{note}
Wave's security scanning applies to any container used in your pipeline, whether it was built by Wave or simply accessed through it. The security scan automatically expires after one week. If a container is accessed again after 7 days or more, the scan will be re-executed.
:::

### Run pipelines using Fusion file system

Wave containers allows you to run your containerised workflow with the {ref}`fusion-page`.

This enables the use of an object storage bucket such as AWS S3 or Google Cloud Storage as your pipeline work directory, simplifying and speeding up many operations on local, AWS Batch, Google Batch or Kubernetes executions.

See the {ref}`Fusion documentation <fusion-page>` for more details.

## Advanced settings

Wave advanced configuration settings are described in the {ref}`Wave <config-wave>` section on the Nextflow configuration page.
