(config-page)=

# Configuration

## Configuration file

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations. Since each configuration file can contain conflicting settings, the sources are ranked to determine which settings are applied. Possible configuration sources, in order of priority:

1. Parameters specified on the command line (`--something value`)
2. Parameters provided using the `-params-file` option
3. Config file specified using the `-c my_config` option
4. The config file named `nextflow.config` in the current directory
5. The config file named `nextflow.config` in the workflow project directory
6. The config file `$HOME/.nextflow/config`
7. Values defined within the pipeline script itself (e.g. `main.nf`)

When more than one of these options for specifying configurations are used, they are merged, so that the settings in the first override the same settings appearing in the second, and so on.

:::{tip}
If you want to ignore any default configuration files and use only a custom one, use `-C <config file>`.
:::

### Config syntax

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax:

```groovy
name = value
```

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (`true`, `false`) do not. Also note that values are typed. This means that, for example, `1` is different from `'1'` — the former is interpreted as the number one, while the latter is interpreted as a string value.

### Config variables

Configuration properties can be used as variables in the configuration file by using the usual `$propertyName` or `${expression}` syntax.

For example:

```groovy
propertyOne = 'world'
anotherProp = "Hello $propertyOne"
customPath = "$PATH:/my/app/folder"
```

Please note, the usual rules for {ref}`string-interpolation` are applied, thus a string containing a variable reference must be wrapped in double-quote chars instead of single-quote chars.

The same mechanism allows you to access environment variables defined in the hosting system. Any variable name not defined in the Nextflow configuration file(s) is interpreted to be a reference to an environment variable with that name. So, in the above example, the property `customPath` is defined as the current system `PATH` to which the string `/my/app/folder` is appended.

### Config comments

Configuration files use the same conventions for comments used by the Groovy or Java programming languages. Thus, use `//` to comment a single line, or `/*` .. `*/` to comment a block on multiple lines.

### Config include

A configuration file can include one or more configuration files using the keyword `includeConfig`. For example:

```groovy
process.executor = 'sge'
process.queue = 'long'
process.memory = '10G'

includeConfig 'path/foo.config'
```

When a relative path is used, it is resolved against the actual location of the including file.

## Config scopes

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope identifier, or grouping the properties in the same scope using the curly brackets notation. For example:

```groovy
alpha.x = 1
alpha.y = 'string value..'

beta {
     p = 2
     q = 'another string ..'
}
```

(config-apptainer)=

### Scope `apptainer`

The `apptainer` scope controls how [Apptainer](https://apptainer.org) containers are executed by Nextflow.

The following settings are available:

`apptainer.autoMounts`
: When `true` Nextflow automatically mounts host paths in the executed container. It requires the `user bind control` feature to be enabled in your Apptainer installation (default: `false`).

`apptainer.cacheDir`
: The directory where remote Apptainer images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`apptainer.enabled`
: Enable Apptainer execution (default: `false`).

`apptainer.engineOptions`
: This attribute can be used to provide any option supported by the Apptainer engine i.e. `apptainer [OPTIONS]`.

`apptainer.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`apptainer.noHttps`
: Pull the Apptainer image with http protocol (default: `false`).

`apptainer.pullTimeout`
: The amount of time the Apptainer pull can last, exceeding which the process is terminated (default: `20 min`).

`apptainer.registry`
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`apptainer.runOptions`
: This attribute can be used to provide any extra command line options supported by `apptainer exec`.

Read the {ref}`container-apptainer` page to learn more about how to use Apptainer containers with Nextflow.

(config-aws)=

### Scope `aws`

The `aws` scope controls the interactions with AWS, including AWS Batch and S3. For example:

```groovy
aws {
    accessKey = '<YOUR S3 ACCESS KEY>'
    secretKey = '<YOUR S3 SECRET KEY>'
    region = 'us-east-1'

    client {
        maxConnections = 20
        connectionTimeout = 10000
        uploadStorageClass = 'INTELLIGENT_TIERING'
        storageEncryption = 'AES256'
    }
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
        maxTransferAttempts = 3
        delayBetweenAttempts = '5 sec'
    }
}
```

:::{tip}
This scope can also be used to configure access to S3-compatible storage outside of AWS, such as [Ceph](https://ceph.com/en/) and [MinIO](https://min.io/).
:::

Read the {ref}`aws-page` and {ref}`amazons3-page` pages for more information.

The following settings are available:

`aws.accessKey`
: AWS account access key

`aws.profile`
: :::{versionadded} 22.12.0-edge
  :::
: AWS profile from `~/.aws/credentials`

`aws.region`
: AWS region (e.g. `us-east-1`)

`aws.secretKey`
: AWS account secret key

`aws.batch.cliPath`
: The path where the AWS command line tool is installed in the host AMI.

`aws.batch.delayBetweenAttempts`
: Delay between download attempts from S3 (default: `10 sec`).

`aws.batch.jobRole`
: The AWS Job Role ARN that needs to be used to execute the Batch Job.

`aws.batch.logsGroup`
: :::{versionadded} 22.09.0-edge
  :::
: The name of the logs group used by Batch Jobs (default: `/aws/batch`).

`aws.batch.maxParallelTransfers`
: Max parallel upload/download transfer operations *per job* (default: `4`).

`aws.batch.maxSpotAttempts`
: :::{versionadded} 22.04.0
  :::
: Max number of execution attempts of a job interrupted by a EC2 spot reclaim event (default: `5`)

`aws.batch.maxTransferAttempts`
: Max number of downloads attempts from S3 (default: `1`).

`aws.batch.retryMode`
: The retry mode configuration setting, to accommodate rate-limiting on [AWS services](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html) (default: `standard`)

`aws.batch.schedulingPriority`
: :::{versionadded} 23.01.0-edge
  :::
: The scheduling priority for all tasks when using [fair-share scheduling for AWS Batch](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/) (default: `0`)

`aws.batch.shareIdentifier`
: :::{versionadded} 22.09.0-edge
  :::
: The share identifier for all tasks when using [fair-share scheduling for AWS Batch](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/)

`aws.batch.volumes`
: One or more container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. `/host/path:/mount/path[:ro|rw]`. Multiple mounts can be specified separating them with a comma or using a list object.

`aws.client.anonymous`
: Allow the access of public S3 buckets without the need to provide AWS credentials. Any service that does not accept unsigned requests will return a service access error.

`aws.client.s3Acl`
: Allow the setting of predefined bucket permissions, also known as *canned ACL*. Permitted values are `Private`, `PublicRead`, `PublicReadWrite`, `AuthenticatedRead`, `LogDeliveryWrite`, `BucketOwnerRead`, `BucketOwnerFullControl`, and `AwsExecRead`. See [Amazon docs](https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl) for details.

`aws.client.connectionTimeout`
: The amount of time to wait (in milliseconds) when initially establishing a connection before timing out.

`aws.client.endpoint`
: The AWS S3 API entry point e.g. `s3-us-west-1.amazonaws.com`.

`aws.client.glacierAutoRetrieval`
: :::{versionadded} 22.12.0-edge
  :::
: *Experimental: may change in a future release.*
: Enable auto retrieval of S3 objects stored with Glacier class store (default: `false`).

`aws.client.glacierExpirationDays`
: :::{versionadded} 22.12.0-edge
  :::
: *Experimental: may change in a future release.*
: The time, in days, between when an object is restored to the bucket and when it expires (default: `7`).

`aws.client.glacierRetrievalTier`
: :::{versionadded} 23.03.0-edge
  :::
: *Experimental: may change in a future release.*
: The retrieval tier to use when restoring objects from Glacier, one of [`Expedited`, `Standard`, `Bulk`].

`aws.client.maxConnections`
: The maximum number of allowed open HTTP connections.

`aws.client.maxErrorRetry`
: The maximum number of retry attempts for failed retryable requests.

`aws.client.protocol`
: The protocol (i.e. HTTP or HTTPS) to use when connecting to AWS.

`aws.client.proxyHost`
: The proxy host to connect through.

`aws.client.proxyPort`
: The port on the proxy host to connect through.

`aws.client.proxyUsername`
: The user name to use when connecting through a proxy.

`aws.client.proxyPassword`
: The password to use when connecting through a proxy.

`aws.client.s3PathStyleAccess`
: Enable the use of path-based access model that is used to specify the address of an object in S3-compatible storage systems.

`aws.client.signerOverride`
: The name of the signature algorithm to use for signing requests made by the client.

`aws.client.socketSendBufferSizeHint`
: The Size hint (in bytes) for the low level TCP send buffer.

`aws.client.socketRecvBufferSizeHint`
: The Size hint (in bytes) for the low level TCP receive buffer.

`aws.client.socketTimeout`
: The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out.

`aws.client.storageEncryption`
: The S3 server side encryption to be used when saving objects on S3, either `AES256` or `aws:kms` values are allowed.

`aws.client.storageKmsKeyId`
: :::{versionadded} 22.05.0-edge
  :::
: The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket ().

`aws.client.userAgent`
: The HTTP user agent header passed with all HTTP requests.

`aws.client.uploadChunkSize`
: The size of a single part in a multipart upload (default: `100 MB`).

`aws.client.uploadMaxAttempts`
: The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).

`aws.client.uploadMaxThreads`
: The maximum number of threads used for multipart upload.

`aws.client.uploadRetrySleep`
: The time to wait after a failed upload attempt to retry the part upload (default: `500ms`).

`aws.client.uploadStorageClass`
: The S3 storage class applied to stored objects, one of \[`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`\] (default: `STANDARD`).

(config-azure)=

### Scope `azure`

The `azure` scope allows you to configure the interactions with Azure, including Azure Batch and Azure Blob Storage.

Read the {ref}`azure-page` page for more information.

The following settings are available:

`azure.activeDirectory.servicePrincipalId`
: The service principal client ID

`azure.activeDirectory.servicePrincipalSecret`
: The service principal client secret

`azure.activeDirectory.tenantId`
: The Azure tenant ID

`azure.batch.accountName`
: The batch service account name.

`azure.batch.accountKey`
: The batch service account key.

`azure.batch.allowPoolCreation`
: Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: `false`).

`azure.batch.autoPoolMode`
: Enable the automatic creation of batch pools depending on the pipeline resources demand (default: `true`).

`azure.batch.copyToolInstallMode`
: Specify where the `azcopy` tool used by Nextflow. When `node` is specified it's copied once during the pool creation. When `task` is provider, it's installed for each task execution (default: `node`).

`azure.batch.terminateJobsOnCompletion`
: Enables the Batch Job to automatically terminate a job once all tasks have completed (default: `true`).

`azure.batch.deleteJobsOnCompletion`
: Enable the automatic deletion of jobs created by the pipeline execution (default: `true`).

`azure.batch.deletePoolsOnCompletion`
: Enable the automatic deletion of compute node pools upon pipeline completion (default: `false`).

`azure.batch.endpoint`
: The batch service endpoint e.g. `https://nfbatch1.westeurope.batch.azure.com`.

`azure.batch.location`
: The name of the batch service region, e.g. `westeurope` or `eastus2`. This is not needed when the endpoint is specified.

`azure.batch.pools.<name>.autoScale`
: Enable autoscaling feature for the pool identified with `<name>`.

`azure.batch.pools.<name>.fileShareRootPath`
: *New in `nf-azure` version `0.11.0`*
: If mounting File Shares, this is the internal root mounting point. Must be `/mnt/resource/batch/tasks/fsmounts` for CentOS nodes or `/mnt/batch/tasks/fsmounts` for Ubuntu nodes (default is for CentOS).

`azure.batch.pools.<name>.maxVmCount`
: Specify the max of virtual machine when using auto scale option.

`azure.batch.pools.<name>.mountOptions`
: *New in `nf-azure` version `0.11.0`*
: Specify the mount options for mounting the file shares (default: `-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp`).

`azure.batch.pools.<name>.offer`
: *New in `nf-azure` version `0.11.0`*
: Specify the offer type of the virtual machine type used by the pool identified with `<name>` (default: `centos-container`).

`azure.batch.pools.<name>.privileged`
: Enable the task to run with elevated access. Ignored if `runAs` is set (default: `false`).

`azure.batch.pools.<name>.publisher`
: *New in `nf-azure` version `0.11.0`*
: Specify the publisher of virtual machine type used by the pool identified with `<name>` (default: `microsoft-azure-batch`).

`azure.batch.pools.<name>.runAs`
: Specify the username under which the task is run. The user must already exist on each node of the pool.

`azure.batch.pools.<name>.scaleFormula`
: Specify the scale formula for the pool identified with `<name>`. See Azure Batch [scaling documentation](https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling) for details.

`azure.batch.pools.<name>.scaleInterval`
: Specify the interval at which to automatically adjust the Pool size according to the autoscale formula. The minimum and maximum value are 5 minutes and 168 hours respectively (default: `10 mins`).

`azure.batch.pools.<name>.schedulePolicy`
: Specify the scheduling policy for the pool identified with `<name>`. It can be either `spread` or `pack` (default: `spread`).

`azure.batch.pools.<name>.sku`
: *New in `nf-azure` version `0.11.0`*
: Specify the ID of the Compute Node agent SKU which the pool identified with `<name>` supports (default: `batch.node.centos 8`).

`azure.batch.pools.<name>.virtualNetwork`
: :::{versionadded} 23.03.0-edge
  :::
: Specify the subnet ID of a virtual network in which to create the pool.

`azure.batch.pools.<name>.vmCount`
: Specify the number of virtual machines provisioned by the pool identified with `<name>`.

`azure.batch.pools.<name>.vmType`
: Specify the virtual machine type used by the pool identified with `<name>`.

`azure.registry.server`
: *New in `nf-azure` version `0.9.8`*
: Specify the container registry from which to pull the Docker images (default: `docker.io`).

`azure.registry.userName`
: *New in `nf-azure` version `0.9.8`*
: Specify the username to connect to a private container registry.

`azure.registry.password`
: *New in `nf-azure` version `0.9.8`*
: Specify the password to connect to a private container registry.

`azure.retryPolicy.delay`
: Delay when retrying failed API requests (default: `500ms`).

`azure.retryPolicy.jitter`
: Jitter value when retrying failed API requests (default: `0.25`).

`azure.retryPolicy.maxAttempts`
: Max attempts when retrying failed API requests (default: `10`).

`azure.retryPolicy.maxDelay`
: Max delay when retrying failed API requests (default: `90s`).

`azure.storage.accountName`
: The blob storage account name

`azure.storage.accountKey`
: The blob storage account key

`azure.storage.sasToken`
: The blob storage shared access signature token. This can be provided as an alternative to the `accountKey` setting.

`azure.storage.tokenDuration`
: The duration of the shared access signature token created by Nextflow when the `sasToken` option is *not* specified (default: `48h`).

(config-charliecloud)=

### Scope `charliecloud`

The `charliecloud` scope controls how [Charliecloud](https://hpc.github.io/charliecloud/) containers are executed by Nextflow.

The following settings are available:

`charliecloud.cacheDir`
: The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`charliecloud.enabled`
: Enable Charliecloud execution (default: `false`).

`charliecloud.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`charliecloud.pullTimeout`
: The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: `20 min`).

`charliecloud.runOptions`
: This attribute can be used to provide any extra command line options supported by the `ch-run` command.

`charliecloud.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

Read the {ref}`container-charliecloud` page to learn more about how to use Charliecloud containers with Nextflow.

(config-conda)=

### Scope `conda`

The `conda` scope controls the creation of a Conda environment by the Conda package manager.

The following settings are available:

`conda.cacheDir`
: Defines the path where Conda environments are stored. When using a compute cluster make sure to provide a shared file system path accessible from all compute nodes.

`conda.createOptions`
: Defines any extra command line options supported by the `conda create` command. For details see the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/commands/create.html).

`conda.createTimeout`
: Defines the amount of time the Conda environment creation can last. The creation process is terminated when the timeout is exceeded (default: `20 min`).

`conda.useMamba`
: Uses the `mamba` binary instead of `conda` to create the Conda environments. For details see the [Mamba documentation](https://github.com/mamba-org/mamba).

`conda.useMicromamba`
: :::{versionadded} 22.05.0-edge
  :::
: uses the `micromamba` binary instead of `conda` to create the Conda environments. For details see the [Micromamba documentation](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

Read the {ref}`conda-page` page to learn more about how to use Conda environments with Nextflow.

(config-dag)=

### Scope `dag`

The `dag` scope controls the layout of the execution graph diagram generated by Nextflow.

The following settings are available:

`dag.enabled`
: When `true` turns on the generation of the DAG file (default: `false`).

`dag.file`
: Graph file name (default: `dag-<timestamp>.dot`).

`dag.overwrite`
: When `true` overwrites any existing DAG file with the same name.

Read the {ref}`dag-visualisation` page to learn more about the execution graph that can be generated by Nextflow.

(config-docker)=

### Scope `docker`

The `docker` scope controls how [Docker](https://www.docker.com) containers are executed by Nextflow.

The following settings are available:

`docker.enabled`
: Enable Docker execution (default: `false`).

`docker.engineOptions`
: This attribute can be used to provide any option supported by the Docker engine i.e. `docker [OPTIONS]`.

`docker.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`docker.fixOwnership`
: Fix ownership of files created by the docker container.

`docker.legacy`
: Use command line options removed since Docker 1.10.0 (default: `false`).

`docker.mountFlags`
: Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`.

`docker.registry`
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`docker.remove`
: Clean-up the container after the execution (default: `true`). See the [Docker documentation](https://docs.docker.com/engine/reference/run/#clean-up---rm) for details.

`docker.runOptions`
: This attribute can be used to provide any extra command line options supported by the `docker run` command. See the [Docker documentation](https://docs.docker.com/engine/reference/run/) for details.

`docker.sudo`
: Executes Docker run command as `sudo` (default: `false`).

`docker.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

`docker.tty`
: Allocates a pseudo-tty (default: `false`).

Read the {ref}`container-docker` page to learn more about how to use Docker containers with Nextflow.

(config-env)=

### Scope `env`

The `env` scope allows the definition one or more variables that will be exported into the environment where workflow tasks are executed.

Simply prefix your variable names with the `env` scope or surround them by curly brackets, as shown below:

```groovy
env.ALPHA = 'some value'
env.BETA = "$HOME/some/path"

env {
     DELTA = 'one more'
     GAMMA = "/my/path:$PATH"
}
```

:::{note}
In the above example, variables like `$HOME` and `$PATH` are evaluated when the workflow is launched. If you want these variables to be evaluated during task execution, escape them with `\$`. This difference is important for variables like `$PATH`, which may be different in the workflow environment versus the task environment.
:::

:::{warning}
The `env` scope provides environment variables to *tasks*, not Nextflow itself. Nextflow environment variables such as `NXF_VER` should be set in the environment in which Nextflow is launched.
:::

(config-executor)=

### Scope `executor`

The `executor` scope controls various executor behaviors.

The following settings are available:

`executor.cpus`
: The maximum number of CPUs made available by the underlying system. Used only by the `local` executor.

`executor.dumpInterval`
: Determines how often to log the executor status (default: `5min`).

`executor.exitReadTimeout`
: Determines how long to wait before returning an error status when a process is terminated but the `.exitcode` file does not exist or is empty (default: `270 sec`). Used only by grid executors.

`executor.jobName`
: Determines the name of jobs submitted to the underlying cluster executor e.g. `executor.jobName = { "$task.name - $task.hash" }`. Make sure the resulting job name matches the validation constraints of the underlying batch scheduler.

`executor.killBatchSize`
: Determines the number of jobs that can be killed in a single command execution (default: `100`).

`executor.memory`
: The maximum amount of memory made available by the underlying system. Used only by the `local` executor.

`executor.name`
: The name of the executor to be used (default: `local`).

`executor.perCpuMemAllocation`
: :::{versionadded} 23.07.0-edge
  :::
: *Used only by the {ref}`slurm-executor` executor.*
: When `true`, specifies memory allocations for SLURM jobs as `--mem-per-cpu <task.memory / task.cpus>` instead of `--mem <task.memory>`.

`executor.perJobMemLimit`
: Specifies Platform LSF *per-job* memory limit mode. See {ref}`lsf-executor`.

`executor.perTaskReserve`
: Specifies Platform LSF *per-task* memory reserve mode. See {ref}`lsf-executor`.

`executor.pollInterval`
: Determines how often to check for process termination. Default varies for each executor (see below).

`executor.queueGlobalStatus`
: :::{versionadded} 23.01.0-edge
  :::
: Determines how job status is retrieved. When `false` only the queue associated with the job execution is queried. When `true` the job status is queried globally i.e. irrespective of the submission queue (default: `false`).

`executor.queueSize`
: The number of tasks the executor will handle in a parallel manner. Default varies for each executor (see below).

`executor.queueStatInterval`
: Determines how often to fetch the queue status from the scheduler (default: `1min`). Used only by grid executors.

`executor.retry.delay`
: :::{versionadded} 22.03.0-edge
  :::
: Delay when retrying failed job submissions (default: `500ms`). Used only by grid executors.

`executor.retry.jitter`
: :::{versionadded} 22.03.0-edge
  :::
: Jitter value when retrying failed job submissions (default: `0.25`). Used only by grid executors.

`executor.retry.maxAttempt`
: :::{versionadded} 22.03.0-edge
  :::
: Max attempts when retrying failed job submissions (default: `3`). Used only by grid executors.

`executor.retry.maxDelay`
: :::{versionadded} 22.03.0-edge
  :::
: Max delay when retrying failed job submissions (default: `30s`). Used only by grid executors.

`executor.retry.reason`
: :::{versionadded} 22.03.0-edge
  :::
: Regex pattern that when verified cause a failed submit operation to be re-tried (default: `Socket timed out`). Used only by grid executors.

`executor.submitRateLimit`
: Determines the max rate of job submission per time unit, for example `'10sec'` (10 jobs per second) or `'50/2min'` (50 jobs every 2 minutes) (default: unlimited).

Some executor settings have different default values depending on the executor.

| Executor       | `queueSize` | `pollInterval` |
| -------------- | ----------- | -------------- |
| AWS Batch      | `1000`      | `10s`          |
| Azure Batch    | `1000`      | `10s`          |
| Google Batch   | `1000`      | `10s`          |
| Grid Executors | `100`       | `5s`           |
| Kubernetes     | `100`       | `5s`           |
| Local          | N/A         | `100ms`        |

The executor settings can be defined as shown below:

```groovy
executor {
    name = 'sge'
    queueSize = 200
    pollInterval = '30 sec'
}
```

When using two (or more) different executors in your pipeline, you can specify their settings separately by prefixing the executor name with the symbol `$` and using it as special scope identifier. For example:

```groovy
executor {
  $sge {
      queueSize = 100
      pollInterval = '30sec'
  }

  $local {
      cpus = 8
      memory = '32 GB'
  }
}
```

The above configuration example can be rewritten using the dot notation as shown below:

```groovy
executor.$sge.queueSize = 100
executor.$sge.pollInterval = '30sec'
executor.$local.cpus = 8
executor.$local.memory = '32 GB'
```

(config-google)=

### Scope `google`

The `google` scope allows you to configure the interactions with Google Cloud, including Google Cloud Batch, Google Life Sciences, and Google Cloud Storage.

Read the {ref}`google-page` page for more information.

The following settings are available:

`google.enableRequesterPaysBuckets`
: When `true` uses the given Google Cloud project ID as the billing project for storage access. This is required when accessing data from *requester pays enabled* buckets. See [Requester Pays on Google Cloud Storage documentation](https://cloud.google.com/storage/docs/requester-pays) (default: `false`).

`google.httpConnectTimeout`
: :::{versionadded} 23.06.0-edge
  :::
: Defines the HTTP connection timeout for Cloud Storage API requests (default: `'60s'`).

`google.httpReadTimeout`
: :::{versionadded} 23.06.0-edge
  :::
: Defines the HTTP read timeout for Cloud Storage API requests (default: `'60s'`).

`google.location`
: The Google Cloud location where jobs are executed (default: `us-central1`).

`google.project`
: The Google Cloud project ID to use for pipeline execution

`google.region`
: *Available only for Google Life Sciences*
: The Google Cloud region where jobs are executed. Multiple regions can be provided as a comma-separated list. Cannot be used with the `google.zone` option. See the [Google Cloud documentation](https://cloud.google.com/compute/docs/regions-zones/) for a list of available regions and zones.

`google.zone`
: *Available only for Google Life Sciences*
: The Google Cloud zone where jobs are executed. Multiple zones can be provided as a comma-separated list. Cannot be used with the `google.region` option. See the [Google Cloud documentation](https://cloud.google.com/compute/docs/regions-zones/) for a list of available regions and zones.

`google.batch.allowedLocations`
: :::{versionadded} 22.12.0-edge
  :::
: Define the set of allowed locations for VMs to be provisioned. See [Google documentation](https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy) for details (default: no restriction).

`google.batch.bootDiskSize`
: Set the size of the virtual machine boot disk, e.g `50.GB` (default: none).

`google.batch.cpuPlatform`
: Set the minimum CPU Platform, e.g. `'Intel Skylake'`. See [Specifying a minimum CPU Platform for VM instances](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications) (default: none).

`google.batch.network`
: Set network name to attach the VM's network interface to. The value will be prefixed with `global/networks/` unless it contains a `/`, in which case it is assumed to be a fully specified network resource URL. If unspecified, the global default network is used.

`google.batch.serviceAccountEmail`
: Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.

`google.batch.spot`
: When `true` enables the usage of *spot* virtual machines or `false` otherwise (default: `false`).

`google.batch.subnetwork`
: Define the name of the subnetwork to attach the instance to must be specified here, when the specified network is configured for custom subnet creation. The value is prefixed with `regions/subnetworks/` unless it contains a `/`, in which case it is assumed to be a fully specified subnetwork resource URL.

`google.batch.usePrivateAddress`
: When `true` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: `false`).

`google.lifeSciences.bootDiskSize`
: Set the size of the virtual machine boot disk e.g `50.GB` (default: none).

`google.lifeSciences.copyImage`
: The container image run to copy input and output files. It must include the `gsutil` tool (default: `google/cloud-sdk:alpine`).

`google.lifeSciences.cpuPlatform`
: Set the minimum CPU Platform e.g. `'Intel Skylake'`. See [Specifying a minimum CPU Platform for VM instances](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications) (default: none).

`google.lifeSciences.debug`
: When `true` copies the `/google` debug directory in that task bucket directory (default: `false`).

`google.lifeSciences.keepAliveOnFailure`
: :::{versionadded} 21.06.0-edge
  :::
: When `true` and a task complete with an unexpected exit status the associated compute node is kept up for 1 hour. This options implies `sshDaemon=true` (default: `false`).

`google.lifeSciences.network`
: :::{versionadded} 21.03.0-edge
  :::
: Set network name to attach the VM's network interface to. The value will be prefixed with `global/networks/` unless it contains a `/`, in which case it is assumed to be a fully specified network resource URL. If unspecified, the global default network is used.

`google.lifeSciences.preemptible`
: When `true` enables the usage of *preemptible* virtual machines or `false` otherwise (default: `true`).

`google.lifeSciences.serviceAccountEmail`
: :::{versionadded} 20.05.0-edge
  :::
: Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.

`google.lifeSciences.subnetwork`
: :::{versionadded} 21.03.0-edge
  :::
: Define the name of the subnetwork to attach the instance to must be specified here, when the specified network is configured for custom subnet creation. The value is prefixed with `regions/subnetworks/` unless it contains a `/`, in which case it is assumed to be a fully specified subnetwork resource URL.

`google.lifeSciences.sshDaemon`
: When `true` runs SSH daemon in the VM carrying out the job to which it's possible to connect for debugging purposes (default: `false`).

`google.lifeSciences.sshImage`
: The container image used to run the SSH daemon (default: `gcr.io/cloud-genomics-pipelines/tools`).

`google.lifeSciences.usePrivateAddress`
: :::{versionadded} 20.03.0-edge
  :::
: When `true` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: `false`).

`google.storage.delayBetweenAttempts`
: :::{versionadded} 21.06.0-edge
  :::
: Delay between download attempts from Google Storage (default `10 sec`).

`google.storage.downloadMaxComponents`
: :::{versionadded} 21.06.0-edge
  :::
: Defines the value for the option `GSUtil:sliced_object_download_max_components` used by `gsutil` for transfer input and output data (default: `8`).

`google.storage.maxParallelTransfers`
: :::{versionadded} 21.06.0-edge
  :::
: Max parallel upload/download transfer operations *per job* (default: `4`).

`google.storage.maxTransferAttempts`
: :::{versionadded} 21.06.0-edge
  :::
: Max number of downloads attempts from Google Storage (default: `1`).

`google.storage.parallelThreadCount`
: :::{versionadded} 21.06.0-edge
  :::
: Defines the value for the option `GSUtil:parallel_thread_count` used by `gsutil` for transfer input and output data (default: `1`).

(config-k8s)=

### Scope `k8s`

The `k8s` scope controls the deployment and execution of workflow applications in a Kubernetes cluster.

The following settings are available:

`k8s.autoMountHostPaths`
: Automatically mounts host paths in the job pods. Only for development purpose when using a single node cluster (default: `false`).

`k8s.computeResourceType`
: :::{versionadded} 22.05.0-edge
  :::
: Define whether use Kubernetes `Pod` or `Job` resource type to carry out Nextflow tasks (default: `Pod`).

`k8s.context`
: Defines the Kubernetes [configuration context name](https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/) to use.

`k8s.fetchNodeName`
: :::{versionadded} 22.05.0-edge
  :::
: If you trace the hostname, activate this option (default: `false`).

`k8s.httpConnectTimeout`
: :::{versionadded} 22.10.0
  :::
: Defines the Kubernetes client request HTTP connection timeout e.g. `'60s'`.

`k8s.httpReadTimeout`
: :::{versionadded} 22.10.0
  :::
: Defines the Kubernetes client request HTTP connection read timeout e.g. `'60s'`.

`k8s.launchDir`
: Defines the path where the workflow is launched and the user data is stored. This must be a path in a shared K8s persistent volume (default: `<volume-claim-mount-path>/<user-name>`).

`k8s.maxErrorRetry`
: :::{versionadded} 22.09.6-edge
  :::
: Defines the Kubernetes API max request retries (default: 4).

`k8s.namespace`
: Defines the Kubernetes namespace to use (default: `default`).

`k8s.pod`
: Allows the definition of one or more pod configuration options such as environment variables, config maps, secrets, etc. It allows the same settings as the {ref}`process-pod` process directive.

`k8s.projectDir`
: Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: `<volume-claim-mount-path>/projects`).

`k8s.pullPolicy`
: Defines the strategy to be used to pull the container image e.g. `pullPolicy: 'Always'`.

`k8s.runAsUser`
: Defines the user ID to be used to run the containers. Shortcut for the `securityContext` option.

`k8s.securityContext`
: Defines the [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) for all pods.

`k8s.serviceAccount`
: Defines the Kubernetes [service account name](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/) to use.

`k8s.storageClaimName`
: The name of the persistent volume claim where store workflow result data.

`k8s.storageMountPath`
: The path location used to mount the persistent volume claim (default: `/workspace`).

`k8s.storageSubPath`
: The path in the persistent volume to be mounted (default: `/`).

`k8s.workDir`
: Defines the path where the workflow temporary data is stored. This must be a path in a shared K8s persistent volume (default: `<user-dir>/work`).

See the {ref}`k8s-page` page for more details.

(config-mail)=

### Scope `mail`

The `mail` scope controls the mail server used to send email notifications.

The following settings are available:

`mail.debug`
: When `true` enables Java Mail logging for debugging purpose.

`mail.from`
: Default email sender address.

`mail.smtp.host`
: Host name of the mail server.

`mail.smtp.port`
: Port number of the mail server.

`mail.smtp.user`
: User name to connect to the mail server.

`mail.smtp.password`
: User password to connect to the mail server.

`mail.smtp.proxy.host`
: Host name of an HTTP web proxy server that will be used for connections to the mail server.

`mail.smtp.proxy.port`
: Port number for the HTTP web proxy server.

`mail.smtp.*`
: Any SMTP configuration property supported by the [Java Mail API](https://javaee.github.io/javamail/), which Nextflow uses to send emails. See the table of available properties [here](https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties).

For example, the following snippet shows how to configure Nextflow to send emails through the [AWS Simple Email Service](https://aws.amazon.com/ses/):

```groovy
mail {
    smtp.host = 'email-smtp.us-east-1.amazonaws.com'
    smtp.port = 587
    smtp.user = '<Your AWS SES access key>'
    smtp.password = '<Your AWS SES secret key>'
    smtp.auth = true
    smtp.starttls.enable = true
    smtp.starttls.required = true
}
```

:::{note}
Some versions of Java (e.g. Java 11 Corretto) do not default to TLS v1.2, and as a result may have issues with 3rd party integrations that enforce TLS v1.2 (e.g. Azure Active Directory OIDC). This problem can be addressed by setting the following config option:

```groovy
mail {
    smtp.ssl.protocols = 'TLSv1.2'
}
```
:::

(config-manifest)=

### Scope `manifest`

The `manifest` scope allows you to define some meta-data information needed when publishing or running your pipeline.

The following settings are available:

`manifest.author`
: Project author name (use a comma to separate multiple names).

`manifest.defaultBranch`
: Git repository default branch (default: `master`).

`manifest.description`
: Free text describing the workflow project.

`manifest.doi`
: Project related publication DOI identifier.

`manifest.homePage`
: Project home page URL.

`manifest.mainScript`
: Project main script (default: `main.nf`).

`manifest.name`
: Project short name.

`manifest.nextflowVersion`
: Minimum required Nextflow version.

  This setting may be useful to ensure that a specific version is used:

  ```groovy
  manifest.nextflowVersion = '1.2.3'        // exact match
  manifest.nextflowVersion = '1.2+'         // 1.2 or later (excluding 2 and later)
  manifest.nextflowVersion = '>=1.2'        // 1.2 or later
  manifest.nextflowVersion = '>=1.2, <=1.5' // any version in the 1.2 .. 1.5 range
  manifest.nextflowVersion = '!>=1.2'       // with ! prefix, stop execution if current version does not match required version.
  ```

`manifest.recurseSubmodules`
: Pull submodules recursively from the Git repository.

`manifest.version`
: Project version number.

The above options can also be specified in a `manifest` block, for example:

```groovy
manifest {
    homePage = 'http://foo.com'
    description = 'Pipeline does this and that'
    mainScript = 'foo.nf'
    version = '1.0.0'
}
```

Read the {ref}`sharing-page` page to learn how to publish your pipeline to GitHub, BitBucket or GitLab.

(config-notification)=

### Scope `notification`

The `notification` scope allows you to define the automatic sending of a notification email message when the workflow execution terminates.

`notification.binding`
: An associative array modelling the variables in the template file.

`notification.enabled`
: Enables the sending of a notification message when the workflow execution completes.

`notification.from`
: Sender address for the notification email message.

`notification.template`
: Path of a template file which provides the content of the notification message.

`notification.to`
: Recipient address for the notification email. Multiple addresses can be specified separating them with a comma.

The notification message is sent my using the STMP server defined in the configuration {ref}`mail scope<config-mail>`.

If no mail configuration is provided, it tries to send the notification message by using the external mail command eventually provided by the underlying system (e.g. `sendmail` or `mail`).

(config-params)=

### Scope `params`

The `params` scope allows you to define parameters that will be accessible in the pipeline script. Simply prefix the parameter names with the `params` scope or surround them by curly brackets, as shown below:

```groovy
params.custom_param = 123
params.another_param = 'string value .. '

params {
    alpha_1 = true
    beta_2 = 'another string ..'
}
```

(config-podman)=

### Scope `podman`

The `podman` scope controls how [Podman](https://podman.io/) containers are executed by Nextflow.

The following settings are available:

`podman.enabled`
: Enable Podman execution (default: `false`).

`podman.engineOptions`
: This attribute can be used to provide any option supported by the Podman engine i.e. `podman [OPTIONS]`.

`podman.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`podman.mountFlags`
: Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`.

`podman.registry`
: The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`podman.remove`
: Clean-up the container after the execution (default: `true`).

`podman.runOptions`
: This attribute can be used to provide any extra command line options supported by the `podman run` command.

`podman.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

Read the {ref}`container-podman` page to learn more about how to use Podman containers with Nextflow.

(config-process)=

### Scope `process`

The `process` scope allows you to specify default {ref}`directives <process-directives>` for processes in your pipeline.

For example:

```groovy
process {
    executor = 'sge'
    queue = 'long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
}
```

By using this configuration, all processes in your pipeline will be executed through the SGE cluster, with the specified settings.

(config-process-selectors)=

#### Process selectors

The `withLabel` selectors allow the configuration of all processes annotated with a {ref}`process-label` directive as shown below:

```groovy
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
        queue = 'long'
    }
}
```

The above configuration example assigns 16 cpus, 64 Gb of memory and the `long` queue to all processes annotated with the `big_mem` label.

In the same manner, the `withName` selector allows the configuration of a specific process in your pipeline by its name. For example:

```groovy
process {
    withName: hello {
        cpus = 4
        memory = 8.GB
        queue = 'short'
    }
}
```

:::{tip}
Label and process names do not need to be enclosed with quotes, provided the name does not include special characters (`-`, `!`, etc) and is not a keyword or a built-in type identifier. When in doubt, you can enclose the label name or process name with single or double quotes.
:::

(config-selector-expressions)=

#### Selector expressions

Both label and process name selectors allow the use of a regular expression in order to apply the same configuration to all processes matching the specified pattern condition. For example:

```groovy
process {
    withLabel: 'foo|bar' {
        cpus = 2
        memory = 4.GB
    }
}
```

The above configuration snippet sets 2 cpus and 4 GB of memory to the processes annotated with a label `foo` and `bar`.

A process selector can be negated prefixing it with the special character `!`. For example:

```groovy
process {
    withLabel: 'foo' { cpus = 2 }
    withLabel: '!foo' { cpus = 4 }
    withName: '!align.*' { queue = 'long' }
}
```

The above configuration snippet sets 2 cpus for the processes annotated with the `foo` label and 4 cpus to all processes *not* annotated with that label. Finally it sets the use of `long` queue to all process whose name does *not* start with `align`.

(config-selector-priority)=

#### Selector priority

When mixing generic process configuration and selectors the following priority rules are applied (from lower to higher):

1. Process generic configuration.
2. Process specific directive defined in the workflow script.
3. `withLabel` selector definition.
4. `withName` selector definition.

For example:

```groovy
process {
    cpus = 4
    withLabel: foo { cpus = 8 }
    withName: bar { cpus = 32 }
}
```

Using the above configuration snippet, all workflow processes use 4 cpus if not otherwise specified in the workflow script. Moreover processes annotated with the `foo` label use 8 cpus. Finally the process named `bar` uses 32 cpus.

(config-report)=

### Scope `report`

The `report` scope allows you to configure the workflow {ref}`execution-report`.

The following settings are available:

`report.enabled`
: If `true` it create the workflow execution report.

`report.file`
: The path of the created execution report file (default: `report-<timestamp>.html`).

`report.overwrite`
: When `true` overwrites any existing report file with the same name.

(config-sarus)=

### Scope `sarus`

The ``sarus`` scope controls how [Sarus](https://sarus.readthedocs.io) containers are executed by Nextflow.

The following settings are available:

`sarus.enabled`
: Enable Sarus execution (default: `false`).

`sarus.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`sarus.runOptions`
: This attribute can be used to provide any extra command line options supported by the `sarus run` command. For details see the [Sarus user guide](https://sarus.readthedocs.io/en/stable/user/user_guide.html).

`sarus.tty`
: Allocates a pseudo-tty (default: `false`).

Read the {ref}`container-sarus` page to learn more about how to use Sarus containers with Nextflow.

(config-shifter)=

### Scope `shifter`

The `shifter` scope controls how [Shifter](https://docs.nersc.gov/programming/shifter/overview/) containers are executed by Nextflow.

The following settings are available:

`shifter.enabled`
: Enable Shifter execution (default: `false`).

Read the {ref}`container-shifter` page to learn more about how to use Shifter containers with Nextflow.

(config-singularity)=

### Scope `singularity`

The `singularity` scope controls how [Singularity](https://sylabs.io/singularity/) containers are executed by Nextflow.

The following settings are available:

`singularity.autoMounts`
: When `true` Nextflow automatically mounts host paths in the executed container. It requires the `user bind control` feature to be enabled in your Singularity installation (default: `false`).

`singularity.cacheDir`
: The directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`singularity.enabled`
: Enable Singularity execution (default: `false`).

`singularity.engineOptions`
: This attribute can be used to provide any option supported by the Singularity engine i.e. `singularity [OPTIONS]`.

`singularity.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`singularity.noHttps`
: Pull the Singularity image with http protocol (default: `false`).

`singularity.pullTimeout`
: The amount of time the Singularity pull can last, exceeding which the process is terminated (default: `20 min`).

`singularity.registry`
: :::{versionadded} 22.12.0-edge
  :::
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`singularity.runOptions`
: This attribute can be used to provide any extra command line options supported by `singularity exec`.

Read the {ref}`container-singularity` page to learn more about how to use Singularity containers with Nextflow.

(config-spack)=

### Scope `spack`

The `spack` scope controls the creation of a Spack environment by the Spack package manager.

The following settings are available:

`spack.cacheDir`
: Defines the path where Spack environments are stored. When using a compute cluster make sure to provide a shared file system path accessible from all compute nodes.

`spack.checksum`
: Enables checksum verification for source tarballs (default: `true`). Only disable when requesting a package version not yet encoded in the corresponding Spack recipe.

`spack.createTimeout`
: Defines the amount of time the Spack environment creation can last. The creation process is terminated when the timeout is exceeded (default: `60 min`).

`spack.parallelBuilds`
: Sets number of parallel package builds (Spack default: coincides with number of available CPU cores).

Nextflow does not allow for fine-grained configuration of the Spack package manager. Instead, this has to be performed directly on the host Spack installation. For more information see the [Spack documentation](https://spack.readthedocs.io).

(config-timeline)=

### Scope `timeline`

The `timeline` scope controls the execution timeline report generated by Nextflow.

The following settings are available:

`timeline.enabled`
: When `true` enables the generation of the timeline report file (default: `false`).

`timeline.file`
: Timeline file name (default: `timeline-<timestamp>.html`).

`timeline.overwrite`
: When `true` overwrites any existing timeline file with the same name.

(config-tower)=

### Scope `tower`

The `tower` scope controls the settings for the [Nextflow Tower](https://tower.nf) monitoring and tracing service.

The following settings are available:

`tower.accessToken`
: The unique access token specific to your account on an instance of Tower.

  Your `accessToken` can be obtained from your Tower instance in the [Tokens page](https://tower.nf/tokens).

`tower.enabled`
: When `true` Nextflow sends the workflow tracing and execution metrics to the Nextflow Tower service (default: `false`).

`tower.endpoint`
: The endpoint of your Tower deployment (default: `https://tower.nf`).

`tower.workspaceId`
: The ID of the Tower workspace where the run should be added (default: the launching user personal workspace).

  The Tower workspace ID can also be specified using the environment variable `TOWER_WORKSPACE_ID` (config file has priority over the environment variable).

(config-trace)=

### Scope `trace`

The `trace` scope controls the layout of the execution trace file generated by Nextflow.

The following settings are available:

`trace.enabled`
: When `true` turns on the generation of the execution trace report file (default: `false`).

`trace.fields`
: Comma separated list of fields to be included in the report. The available fields are listed at {ref}`this page <trace-fields>`.

`trace.file`
: Trace file name (default: `trace-<timestamp>.txt`).

`trace.overwrite`
: When `true` overwrites any existing trace file with the same name.

`trace.raw`
: When `true` turns on raw number report generation i.e. date and time are reported as milliseconds and memory as number of bytes.

`trace.sep`
: Character used to separate values in each row (default: `\t`).

The above options can also be specified in a `trace` block, for example:

```groovy
trace {
    enabled = true
    file = 'pipeline_trace.txt'
    fields = 'task_id,name,status,exit,realtime,%cpu,rss'
}
```

Read the {ref}`trace-report` page to learn more about the execution report that can be generated by Nextflow.

(config-miscellaneous)=

### Miscellaneous

There are additional variables that can be defined within a configuration file that do not have a dedicated scope.

`cleanup`
: If `true`, on a successful completion of a run all files in *work* directory are automatically deleted.

  :::{warning}
  The use of the `cleanup` option will prevent the use of the *resume* feature on subsequent executions of that pipeline run. Also, be aware that deleting all scratch files can take a lot of time, especially when using a shared file system or remote cloud storage.
  :::

`dumpHashes`
: If `true`, dump task hash keys in the log file, for debugging purposes.

(config-profiles)=

## Config profiles

Configuration files can contain the definition of one or more *profiles*. A profile is a set of configuration attributes that can be selected during pipeline execution by using the `-profile` command line option.

Configuration profiles are defined by using the special scope `profiles`, which group the attributes that belong to the same profile using a common prefix. For example:

```groovy
profiles {

    standard {
        process.executor = 'local'
    }

    cluster {
        process.executor = 'sge'
        process.queue = 'long'
        process.memory = '10GB'
    }

    cloud {
        process.executor = 'cirrus'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }

}
```

This configuration defines three different profiles: `standard`, `cluster`, and `cloud`, that each set different process
configuration strategies depending on the target runtime platform. The `standard` profile is used by default when no profile is specified.

:::{tip}
Multiple configuration profiles can be specified by separating the profile names with a comma, for example:

```bash
nextflow run <your script> -profile standard,cloud
```
:::

:::{danger}
When using the `profiles` feature in your config file, do NOT set attributes in the same scope both inside and outside a `profiles` context. For example:

```groovy
process.cpus = 1

profiles {
  foo {
    process.memory = '2 GB'
  }

  bar {
    process.memory = '4 GB'
  }
}
```

In the above example, the `process.cpus` attribute is not correctly applied because the `process` scope is also used in the `foo` and `bar` profiles.
:::

(config-env-vars)=

## Environment variables

The following environment variables control the configuration of the Nextflow runtime and the underlying Java virtual machine.

`NXF_ANSI_LOG`
: Enables/disables ANSI console output (default `true` when ANSI terminal is detected).

`NXF_ANSI_SUMMARY`
: Enables/disables ANSI completion summary: `true\|false` (default: print summary if execution last more than 1 minute).

`NXF_ASSETS`
: Defines the directory where downloaded pipeline repositories are stored (default: `$NXF_HOME/assets`)

`NXF_CLOUDCACHE_PATH`
: :::{versionadded} 23.07.0-edge
  :::
: Defines the base cache path when using the cloud cache store.

`NXF_CHARLIECLOUD_CACHEDIR`
: Directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CLASSPATH`
: Allows the extension of the Java runtime classpath with extra JAR files or class folders.

`NXF_CLOUD_DRIVER`
: Defines the default cloud driver to be used if not specified in the config file or as command line option, either `aws` or `google`.

`NXF_CONDA_CACHEDIR`
: Directory where Conda environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CONDA_ENABLED`
: :::{versionadded} 22.08.0-edge
  :::
: Enable the use of Conda recipes defined by using the {ref}`process-conda` directive. (default: `false`).

`NXF_DEBUG`
: Defines scripts debugging level: `1` dump task environment variables in the task log file; `2` enables command script execution tracing; `3` enables command wrapper execution tracing.

`NXF_DEFAULT_DSL`
: :::{versionadded} 22.03.0-edge
  :::
: Defines the DSL version that should be used in not specified otherwise in the script of config file (default: `2`)

`NXF_DISABLE_JOBS_CANCELLATION`
: :::{versionadded} 21.12.0-edge
  :::
: Disables the cancellation of child jobs on workflow execution termination.

`NXF_DISABLE_PARAMS_TYPE_DETECTION`
: :::{versionadded} 23.07.0-edge
:::
: Disables the automatic type detection of command line parameters.

`NXF_DISABLE_WAVE_REQUIREMENT`
: :::{versionadded} 23.07.0-edge
:::
: Disables the requirement for Wave service when enabling the Fusion file system.

`NXF_ENABLE_SECRETS`
: :::{versionadded} 21.09.0-edge
  :::
: Enable Nextflow secrets features (default: `true`)

`NXF_ENABLE_STRICT`
: :::{versionadded} 22.05.0-edge
  :::
: Enable Nextflow *strict* execution mode (default: `false`)

`NXF_EXECUTOR`
: Defines the default process executor e.g. `sge`

`NXF_FILE_ROOT`
: :::{versionadded} 23.05.0-edge
  :::
: The file storage path against which relative file paths are resolved.
: For example, with `NXF_FILE_ROOT=/some/root/path`, the use of `file('foo')` will be resolved to the absolute path `/some/root/path/foo`. A remote root path can be specified using the usual protocol prefix, e.g. `NXF_FILE_ROOT=s3://my-bucket/data`. Files defined using an absolute path are not affected by this setting.

`NXF_HOME`
: Nextflow home directory (default: `$HOME/.nextflow`).

`NXF_JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the `JAVA_HOME` variable if defined.

`NXF_JVM_ARGS`
: :::{versionadded} 21.12.1-edge
  :::
: Allows the setting Java VM options. This is similar to `NXF_OPTS` however it's only applied the JVM running Nextflow and not to any java pre-launching commands.

`NXF_OFFLINE`
: When `true` disables the project automatic download and update from remote repositories (default: `false`).

`NXF_OPTS`
: Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of `-Dkey[=value]` properties.

`NXF_ORG`
: Default `organization` prefix when looking for a hosted repository (default: `nextflow-io`).

`NXF_PARAMS_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the pipeline parameters file .

`NXF_PID_FILE`
: Name of the file where the process PID is saved when Nextflow is launched in background.

`NXF_SCM_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the SCM config file .

`NXF_SINGULARITY_CACHEDIR`
: Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SINGULARITY_LIBRARYDIR`
: :::{versionadded} 21.09.0-edge
  :::
: Directory where remote Singularity images are retrieved. It should be a directory accessible to all compute nodes.

`NXF_SPACK_CACHEDIR`
: Directory where Spack environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SPACK_ENABLED`
: :::{versionadded} 23.02.0-edge
  :::
: Enable the use of Spack recipes defined by using the {ref}`process-spack` directive. (default: `false`).

`NXF_TEMP`
: Directory where temporary files are stored

`NXF_VER`
: Defines which version of Nextflow to use.

`NXF_WORK`
: Directory where working files are stored (usually your *scratch* directory)

`NXF_WRAPPER_STAGE_FILE_THRESHOLD`
: :::{versionadded} 23.05.0-edge
  :::
: Defines the minimum size of the `.command.run` staging script for it to be written to a separate `.command.stage` file (default: `'1 MB'`).
: This setting is useful for executors that impose a size limit on job scripts.

`JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow.

`JAVA_CMD`
: Defines the path location of the Java binary command used to launch Nextflow.

`HTTP_PROXY`
: Defines the HTTP proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `http://user:password@proxy-host.com:port`.
  :::

`HTTPS_PROXY`
: Defines the HTTPS proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `https://user:password@proxy-host.com:port`.
  :::

`FTP_PROXY`
: :::{versionadded} 21.06.0-edge
  :::
: Defines the FTP proxy server. Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `ftp://user:password@proxy-host.com:port`.

`NO_PROXY`
: Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.

(config-feature-flags)=

## Feature flags

Some features can be enabled using the `nextflow.enable` and `nextflow.preview` flags. These flags can be specified in the pipeline script or the configuration file, and they are generally used to introduce experimental or other opt-in features.

`nextflow.enable.configProcessNamesValidation`

: When `true`, prints a warning for every `withName:` process selector that doesn't match a process in the pipeline (default: `true`).

`nextflow.enable.dsl`

: Defines the DSL version to use (`1` or `2`).

: :::{versionadded} 22.03.0-edge
  DSL2 is the default DSL version.
  :::

: :::{versionadded} 22.12.0-edge
  DSL1 is no longer supported.
  :::

`nextflow.enable.moduleBinaries`

: When `true`, enables the use of modules with binary scripts. See {ref}`module-binaries` for more information.

`nextflow.enable.strict`

: :::{versionadded} 22.05.0-edge
  :::

: When `true`, the pipeline is executed in "strict" mode, which introduces the following rules:

  - When reading a params file, Nextflow will fail if a dynamic param value references an undefined variable

  - When merging params from a config file with params from the command line, Nextflow will fail if a param is specified from both sources but with different types

  - When using the `join` operator, the `failOnDuplicate` option is `true` by default

  - When using the `join` operator, the `failOnMismatch` option is `true` by default (unless `remainder` is also `true`)

  - When using the `publishDir` process directive, the `failOnError` option is `true` by default

  - In a process definition, Nextflow will fail if an input or output tuple has only one element

  - In a process definition, Nextflow will fail if an output emit name is not a valid identifier (i.e. it should match the pattern `/[A-Za-z_][A-Za-z0-9_]*/`)

  - During a process execution, Nextflow will fail if a received input tuple does not have the same number of elements as was declared

  - During a process execution, Nextflow will fail if the `storeDir` directive is used with non-file outputs

  - Nextflow will fail if a pipeline param is referenced before it is defined

  - Nextflow will fail if multiple functions and/or processes with the same name are defined in a module script

`nextflow.preview.recursion`

: :::{versionadded} 21.11.0-edge
  :::

: *Experimental: may change in a future release.*

: When `true`, enables process and workflow recursion. See [this GitHub discussion](https://github.com/nextflow-io/nextflow/discussions/2521) for more information.
