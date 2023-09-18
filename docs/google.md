(google-page)=

# Google Cloud

## Credentials

Credentials for submitting requests to the Google Cloud Batch and Cloud LifeSciences API are picked up from your environment using [Application Default Credentials](https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http). Application Default Credentials are designed to use the credentials most natural to the environment in which a tool runs.

The most common case will be to pick up your end-user Google credentials from your workstation. You can create these by running the command:

```bash
gcloud auth application-default login
```

and running through the authentication flow. This will write a credential file to your gcloud configuration directory that will be used for any tool you run on your workstation that picks up default credentials.

The next most common case would be when running on a Compute Engine VM. In this case, Application Default Credentials will pick up the Compute Engine Service Account credentials for that VM.

See the [Application Default Credentials](https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http) documentation for how to enable other use cases.

Finally, the `GOOGLE_APPLICATION_CREDENTIALS` environment variable can be used to specify location of the Google credentials file.

If you don't have it, the credentials file can be downloaded from the Google Cloud Console following these steps:

- Open the [Google Cloud Console](https://console.cloud.google.com)
- Go to APIs & Services → Credentials
- Click on the *Create credentials* (blue) drop-down and choose *Service account key*, in the following page
- Select an existing *Service account* or create a new one if needed
- Select JSON as *Key type*
- Click the *Create* button and download the JSON file giving a name of your choice e.g. `creds.json`.

Then, define the following variable replacing the path in the example with the one of your credentials file just downloaded:

```bash
export GOOGLE_APPLICATION_CREDENTIALS="/path/your/file/creds.json"
```

(google-batch)=

## Cloud Batch

:::{versionadded} 22.07.1-edge
:::

[Google Cloud Batch](https://cloud.google.com/batch) is a managed computing service that allows the execution of containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Google Cloud Batch, allowing the seamless deployment of Nextflow pipelines in the cloud, in which tasks are offloaded to the Cloud Batch service.

Read the {ref}`Google Cloud Batch executor <google-batch-executor>` section to learn more about the `google-batch` executor in Nextflow.

(google-batch-config)=

### Configuration

Make sure to have defined in your environment the `GOOGLE_APPLICATION_CREDENTIALS` variable. See the [Credentials](#credentials) section for details.

:::{note}
Make sure your Google account is allowed to access the Google Cloud Batch service by checking the [APIs & Services](https://console.cloud.google.com/apis/dashboard) dashboard.
:::

Create or edit the file `nextflow.config` in your project root directory. The config must specify the following parameters:

- Google Cloud Batch as Nextflow executor
- The Docker container image(s) for pipeline tasks
- The Google Cloud project ID and location

Example:

```groovy
process {
    executor = 'google-batch'
    container = 'your/container:latest'
}

google {
    project = 'your-project-id'
    location = 'us-central1'
}
```

Notes:

- A container image must be specified to execute processes. You can use a different Docker image for each process using one or more {ref}`config-process-selectors`.
- Make sure to specify the project ID, not the project name.
- Make sure to specify a location where Google Batch is available. Refer to the [Google Batch documentation](https://cloud.google.com/batch/docs/get-started#locations) for region availability.

Read the {ref}`Google configuration<config-google>` section to learn more about advanced configuration options.

(google-batch-process)=

### Process definition

Processes can be defined as usual and by default the `cpus` and `memory` directives are used to find the cheapest machine type available at current location that fits the requested resources. If `memory` is not specified, 1GB of memory is allocated per cpu.

:::{versionadded} 23.02.0-edge
The `machineType` directive can be a list of patterns separated by comma. The pattern can contain a `*` to match any number of characters and `?` to match any single character. Examples of valid patterns: `c2-*`, `m?-standard*`, `n*`.

Alternatively it can also be used to define a specific predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) or a custom machine type.
:::

Examples:

```groovy
process automatic_resources_task {
    cpus 8
    memory '40 GB'

    """
    <Your script here>
    """
}

process allowing_some_series {
    cpus 8
    memory '20 GB'
    machineType 'n2-*,c2-*,m3-*'

    """
    <Your script here>
    """
}

process predefined_resources_task {
    machineType 'n1-highmem-8'

    """
    <Your script here>
    """
}
```

:::{versionadded} 23.06.0-edge
:::

The `disk` directive can be used to set the boot disk size or provision a disk for scratch storage. If the disk type is specified with the `type` option, a new disk will be mounted to the task VM at `/tmp` with the requested size and type. Otherwise, it will set the boot disk size, overriding the `google.batch.bootDiskSize` config option. See the [Google Batch documentation](https://cloud.google.com/compute/docs/disks) for more information about the available disk types.

Examples:

```groovy
// set the boot disk size
disk 100.GB

// mount a persistent disk at '/tmp'
disk 100.GB, type: 'pd-standard'

// mount a local SSD disk at '/tmp' (should be a multiple of 375 GB)
disk 375.GB, type: 'local-ssd'
```

### Pipeline execution

The pipeline can be launched either in a local computer or a cloud instance. Pipeline input data can be stored either locally or in a Google Storage bucket.

The pipeline execution must specify a Google Storage bucket where the workflow's intermediate results are stored using the `-work-dir` command line options. For example:

```bash
nextflow run <script or project name> -work-dir gs://my-bucket/some/path
```

:::{tip}
Any input data **not** stored in a Google Storage bucket will automatically be transferred to the pipeline work bucket. Use this feature with caution being careful to avoid unnecessary data transfers.
:::

:::{warning}
The Google Storage path needs to contain at least sub-directory. Don't use only the bucket name e.g. `gs://my-bucket`.
:::

### Spot instances

Spot instances are supported adding the following setting in the Nextflow config file:

```groovy
google {
    batch.spot = true
}
```

Since this type of virtual machines can be retired by the provider before the job completion, it is advisable to add the following retry strategy to your config file to instruct Nextflow to automatically re-execute a job if the virtual machine was terminated preemptively:

```groovy
process {
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
}
```

### Fusion file system

:::{versionadded} 23.02.0-edge
:::

The Google Batch executor supports the use of {ref}`fusion-page`. Fusion allows the use of Google Cloud Storage as a virtual distributed file system, optimising the data transfer and speeding up most job I/O operations.

To enable the use of Fusion file system in your pipeline, add the following snippet to your Nextflow configuration file:

```groovy
fusion.enabled = true
wave.enabled = true
process.scratch = false
tower.accessToken = '<YOUR ACCESS TOKEN>'
```

The [Tower](https://cloud.tower.nf) access token is optional, but it enables higher API rate limits for the {ref}`wave-page` service required by Fusion.

By default, Fusion mounts a local SSD disk to the VM at `/tmp`, using a machine type that can attach local SSD disks. If you specify your own machine type or machine series, they should be able to attach local SSD disks, otherwise the job scheduling will fail.

:::{versionadded} 23.06.0-edge
:::

The `disk` directive can be used to override the disk requested by Fusion. See the {ref}`Process definition <google-batch-process>` section above for examples. Note that local SSD disks must be a multiple of 375 GB in size, otherwise the size will be increased to the next multiple of 375 GB.

### Supported directives

The integration with Google Batch is a developer preview feature. Currently, the following Nextflow directives are supported:

- {ref}`process-accelerator`
- {ref}`process-container`
- {ref}`process-containeroptions`
- {ref}`process-cpus`
- {ref}`process-disk`
- {ref}`process-executor`
- {ref}`process-machinetype`
- {ref}`process-memory`
- {ref}`process-time`

(google-lifesciences)=

## Cloud Life Sciences

:::{versionadded} 20.01.0-edge
:::

:::{note}
In versions of Nextflow prior to `21.04.0`, the following variables must be defined in your system environment:

```bash
export NXF_VER=20.01.0
export NXF_MODE=google
```
:::

[Cloud Life Sciences](https://cloud.google.com/life-sciences/) is a managed computing service that allows the execution of containerized workloads in the Google Cloud Platform infrastructure.

Nextflow provides built-in support for Cloud Life Sciences, allowing the seamless deployment of Nextflow pipelines in the cloud, in which tasks are offloaded to the Cloud Life Sciences service.

Read the {ref}`Google Life Sciences executor <google-lifesciences-executor>` page to learn about the `google-lifesciences` executor in Nextflow.

:::{warning}
This API works well for coarse-grained workloads (i.e. long-running jobs). It's not suggested the use this feature for pipelines spawning many short lived tasks.
:::

(google-lifesciences-config)=

### Configuration

Make sure to have defined in your environment the `GOOGLE_APPLICATION_CREDENTIALS` variable. See the section [Credentials](#credentials) for details.

:::{tip}
Make sure to enable the Cloud Life Sciences API beforehand. To learn how to enable it follow [this link](https://cloud.google.com/life-sciences/docs/quickstart).
:::

Create a `nextflow.config` file in the project root directory. The config must specify the following parameters:

- Google Life Sciences as Nextflow executor
- The Docker container image(s) for pipeline tasks
- The Google Cloud project ID
- The Google Cloud region or zone where the Compute Engine VMs will be executed.
  You need to specify one or the other, *not* both. Multiple regions or zones can be specified as a comma-separated list, e.g. `google.zone = 'us-central1-f,us-central-1-b'`.

Example:

```groovy
process {
    executor = 'google-lifesciences'
    container = 'your/container:latest'
}

google {
    project = 'your-project-id'
    zone = 'europe-west1-b'
}
```

Notes:
- A container image must be specified to execute processes. You can use a different Docker image for each process using one or more {ref}`config-process-selectors`.
- Make sure to specify the project ID, not the project name.
- Make sure to specify a location where Google Life Sciences is available. Refer to the [Google Cloud documentation](https://cloud.google.com/life-sciences/docs/concepts/locations) for details.

Read the {ref}`Google configuration<config-google>` section to learn more about advanced configuration options.

### Process definition

Processes can be defined as usual and by default the `cpus` and `memory` directives are used to instantiate a custom machine type with the specified compute resources. If `memory` is not specified, 1GB of memory is allocated per cpu. A persistent disk will be created with size corresponding to the `disk` directive. If `disk` is not specified, the instance default is chosen to ensure reasonable I/O performance.

The process `machineType` directive may optionally be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) If specified, this value overrides the `cpus` and `memory` directives. If the `cpus` and `memory` directives are used, the values must comply with the allowed custom machine type [specifications](https://cloud.google.com/compute/docs/instances/creating-instance-with-custom-machine-type#specifications) . Extended memory is not directly supported, however high memory or cpu predefined instances may be utilized using the `machineType` directive

Examples:

```groovy
process custom_resources_task {
    cpus 8
    memory '40 GB'
    disk '200 GB'

    """
    <Your script here>
    """
}

process predefined_resources_task {
    machineType 'n1-highmem-8'

    """
    <Your script here>
    """
}
```

### Pipeline execution

The pipeline can be launched either in a local computer or a cloud instance. Pipeline input data can be stored either locally or in a Google Storage bucket.

The pipeline execution must specify a Google Storage bucket where the workflow's intermediate results are stored using the `-work-dir` command line options. For example:

```bash
nextflow run <script or project name> -work-dir gs://my-bucket/some/path
```

:::{tip}
Any input data *not* stored in a Google Storage bucket will be automatically transferred to the pipeline work bucket. Use this feature with caution, being careful to avoid unnecessary data transfers.
:::

### Preemptible instances

Preemptible instances are supported adding the following setting in the Nextflow config file:

```groovy
google {
    lifeSciences.preemptible = true
}
```

Since this type of virtual machines can be retired by the provider before the job completion, it is advisable to add the following retry strategy to your config file to instruct Nextflow to automatically re-execute a job if the virtual machine was terminated preemptively:

```groovy
process {
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
}
```

:::{note}
Preemptible instances have a [runtime limit](https://cloud.google.com/compute/docs/instances/preemptible) of 24 hours.
:::

:::{tip}
For an exhaustive list of error codes, refer to the official Google Life Sciences [documentation](https://cloud.google.com/life-sciences/docs/troubleshooting#error_codes).
:::

### Hybrid execution

Nextflow allows the use of multiple executors in the same workflow. This feature enables the deployment of hybrid workloads, in which some jobs are executed in the local computer or local computing cluster, and some jobs are offloaded to Google Life Sciences.

To enable this feature, use one or more {ref}`config-process-selectors` in your Nextflow configuration file to apply the Google Life Sciences executor to the subset of processes that you want to offload. For example:

```groovy
process {
    withLabel: bigTask {
        executor = 'google-lifesciences'
        container = 'my/image:tag'
    }
}

google {
    project = 'your-project-id'
    zone = 'europe-west1-b'
}
```

Then launch the pipeline with the `-bucket-dir` option to specify a Google Storage path for the jobs computed with Google Life Sciences and, optionally, the `-work-dir` to specify the local storage for the jobs computed locally:

```bash
nextflow run <script or project name> -bucket-dir gs://my-bucket/some/path
```

:::{warning}
The Google Storage path needs to contain at least one sub-directory (e.g. `gs://my-bucket/work` rather than `gs://my-bucket`).
:::

### Limitations

- Compute resources in Google Cloud are subject to [resource quotas](https://cloud.google.com/compute/quotas), which may affect your ability to run pipelines at scale. You can request quota increases, and your quotas may automatically increase over time as you use the platform. In particular, GPU quotas are initially set to 0, so you must explicitly request a quota increase in order to use GPUs. You can initially request an increase to 1 GPU at a time, and after one billing cycle you may be able to increase it further.

- Currently, it's not possible to specify a disk type different from the default one assigned by the service depending on the chosen instance type.

### Troubleshooting

- Make sure to enable the Compute Engine API, Life Sciences API and Cloud Storage API in the [APIs & Services Dashboard](https://console.cloud.google.com/apis/dashboard) page.

- Make sure to have enough compute resources to run your pipeline in your project [Quotas](https://console.cloud.google.com/iam-admin/quotas) (i.e. Compute Engine CPUs, Compute Engine Persistent Disk, Compute Engine In-use IP addresses, etc).

- Make sure your security credentials allow you to access any Google Storage bucket where input data and temporary files are stored.

- When a job fails, you can check the `google/` directory in the task work directory (in the bucket storage), which contains useful information about the job execution. To enable the creation of this directory, set `google.lifeSciences.debug = true` in the Nextflow config.

- You can enable the optional SSH daemon in the job VM by setting `google.lifeSciences.sshDaemon = true` in the Nextflow config.

- Make sure you are choosing a `location` where the [Cloud Life Sciences API is available](https://cloud.google.com/life-sciences/docs/concepts/locations), and a `region` or `zone` where the [Compute Engine API is available](https://cloud.google.com/compute/docs/regions-zones/).
