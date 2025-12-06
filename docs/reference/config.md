(config-options)=

# Configuration options

This page lists all of the available settings in the {ref}`Nextflow configuration <config-page>`.

(config-unscoped)=

## Unscoped options

The following settings are available:

`bucketDir`
: The remote work directory used by hybrid workflows. Equivalent to the `-bucket-dir` option of the `run` command.

`cleanup`
: Delete all files associated with a run in the work directory when the run completes successfully (default: `false`).

: :::{warning}
  This option will prevent the use of the *resume* feature on subsequent executions of that pipeline run.
  :::
  :::{warning}
  This option is not supported for remote work directories, such as Amazon S3, Google Cloud Storage, and Azure Blob Storage.
  :::

`outputDir`
: :::{versionadded} 24.10.0
  :::
: The pipeline output directory. Equivalent to the `-output-dir` option of the `run` command.

`resume`
: Enable the use of previously cached task executions. Equivalent to the `-resume` option of the `run` command.

`workDir`
: The pipeline work directory. Equivalent to the `-work-dir` option of the `run` command.

(config-apptainer)=

## `apptainer`

The `apptainer` scope controls how [Apptainer](https://apptainer.org) containers are executed by Nextflow.

The following settings are available:

`apptainer.autoMounts`
: Automatically mount host paths in the executed container (default: `true`). It requires the `user bind control` feature to be enabled in your Apptainer installation.

`apptainer.cacheDir`
: The directory where remote Apptainer images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`apptainer.enabled`
: Execute tasks with Apptainer containers (default: `false`).

`apptainer.engineOptions`
: Specify additional options supported by the Apptainer engine i.e. `apptainer [OPTIONS]`.

`apptainer.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`apptainer.libraryDir`
: Directory where remote Apptainer images are retrieved. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`apptainer.noHttps`
: Pull the Apptainer image with http protocol (default: `false`).

`apptainer.ociAutoPull`
: :::{versionadded} 23.12.0-edge
  :::
: When enabled, OCI (and Docker) container images are pulled and converted to the SIF format by the Apptainer run command, instead of Nextflow (default: `false`).

  :::{note}
  Leave `ociAutoPull` disabled if you are willing to build a Singularity/Apptainer native image with Wave (see the {ref}`wave-singularity` section).
  :::

`apptainer.pullTimeout`
: The amount of time the Apptainer pull can last, exceeding which the process is terminated (default: `20 min`).

`apptainer.registry`
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`apptainer.runOptions`
: Specify extra command line options supported by `apptainer exec`.

(config-aws)=

## `aws`

The `aws` scope controls the interactions with AWS, including AWS Batch and S3.

The following settings are available:

`aws.accessKey`
: AWS account access key.

`aws.profile`
: :::{versionadded} 22.12.0-edge
  :::
: AWS profile from `~/.aws/credentials`.

`aws.region`
: AWS region (e.g. `us-east-1`).

`aws.secretKey`
: AWS account secret key.

`aws.batch.cliPath`
: The path where the AWS command line tool is installed in the host AMI.

`aws.batch.delayBetweenAttempts`
: Delay between download attempts from S3 (default: `10 sec`).

`aws.batch.executionRole`
: The AWS Batch [Execution Role](https://docs.aws.amazon.com/batch/latest/userguide/execution-IAM-role.html) ARN that needs to be used to execute the Batch Job. It is mandatory when using AWS Fargate.

`aws.batch.jobRole`
: The AWS Batch Job Role ARN that needs to be used to execute the Batch Job.

`aws.batch.logsGroup`
: :::{versionadded} 22.09.0-edge
  :::
: The name of the logs group used by Batch Jobs (default: `/aws/batch/job`).

`aws.batch.maxParallelTransfers`
: Max parallel upload/download transfer operations *per job* (default: `4`).

`aws.batch.maxSpotAttempts`
: :::{versionadded} 22.04.0
  :::
: :::{versionchanged} 24.08.0-edge
  The default value was changed from `5` to `0`.
  :::
: Max number of execution attempts of a job interrupted by a EC2 Spot reclaim event (default: `0`)

`aws.batch.maxTransferAttempts`
: Max number of downloads attempts from S3 (default: `1`).

`aws.batch.platformType`
: The compute platform type used by AWS Batch. Can be either `ec2` or `fargate`. Set to `fargate` to use [AWS Fargate](https://docs.aws.amazon.com/batch/latest/userguide/fargate.html).

`aws.batch.retryMode`
: The [retry mode](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html) used to handle rate-limiting by AWS APIs. Can be one of `standard`, `legacy`, `adaptive`, or `built-in` (default: `standard`).

`aws.batch.schedulingPriority`
: :::{versionadded} 23.01.0-edge
  :::
: The scheduling priority for all tasks when using [fair-share scheduling](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/) (default: `0`).

`aws.batch.shareIdentifier`
: :::{versionadded} 22.09.0-edge
  :::
: The share identifier for all tasks when using [fair-share scheduling](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/).

`aws.batch.terminateUnschedulableJobs`
: :::{versionadded} 25.03.0-edge
  :::
: When `true`, jobs that cannot be scheduled due to lack of resources or misconfiguration are terminated and handled as task failures (default: `false`).

`aws.batch.volumes`
: List of container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. `/host/path:/mount/path[:ro|rw]`.

`aws.client.anonymous`
: Allow the access of public S3 buckets without providing AWS credentials (default: `false`). Any service that does not accept unsigned requests will return a service access error.

`aws.client.connectionTimeout`
: The amount of time to wait (in milliseconds) when initially establishing a connection before timing out (default: `10000`).

`aws.client.endpoint`
: The AWS S3 API entry point e.g. `https://s3-us-west-1.amazonaws.com`. The endpoint must include the protocol prefix e.g. `https://`.

`aws.client.maxConnections`
: The maximum number of open HTTP connections used by the S3 client (default: `50`).

`aws.client.maxDownloadHeapMemory`
: The maximum size for the heap memory buffer used by concurrent downloads. It must be at least 10 times the `minimumPartSize` (default:`400 MB`).

`aws.client.maxErrorRetry`
: The maximum number of retry attempts for failed retryable requests (default: `-1`).

`aws.client.minimumPartSize`
: :::{versionadded} 25.06.0-edge
  :::
: The minimum part size used for multipart S3 transfers (default: `8 MB`).

`aws.client.multipartThreshold`
: :::{versionadded} 25.06.0-edge
  :::
: The object size threshold used for multipart S3 transfers (default: same as `aws.cllient.minimumPartSize`).

`aws.client.protocol`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The protocol to use when connecting to AWS. Can be `http` or `https` (default: `'https'`).

`aws.client.proxyHost`
: The proxy host to connect through.

`aws.client.proxyPassword`
: The password to use when connecting through a proxy.

`aws.client.proxyPort`
: The port to use when connecting through a proxy.

`aws.client.proxyScheme`
: :::{versionadded} 25.06.0-edge
  :::
: The protocol scheme to use when connecting through a proxy. Can be `http` or `https` (default: `'http'`).

`aws.client.proxyUsername`
: The user name to use when connecting through a proxy.

`aws.client.requesterPays`
: :::{versionadded} 24.05.0-edge
  :::
: Use [Requester Pays](https://docs.aws.amazon.com/AmazonS3/latest/userguide/RequesterPaysBuckets.html) for S3 buckets (default: `false`).

`aws.client.s3Acl`
: Specify predefined bucket permissions, also known as [canned ACL](https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl). Can be one of `Private`, `PublicRead`, `PublicReadWrite`, `AuthenticatedRead`, `LogDeliveryWrite`, `BucketOwnerRead`, `BucketOwnerFullControl`, or `AwsExecRead`.

`aws.client.s3PathStyleAccess`
: Use the path-based access model to access objects in S3-compatible storage systems (default: `false`).

`aws.client.signerOverride`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The name of the signature algorithm to use for signing requests made by the client.

`aws.client.socketSendBufferSizeHint`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The Size hint (in bytes) for the low level TCP send buffer (default: `0`).

`aws.client.socketRecvBufferSizeHint`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The Size hint (in bytes) for the low level TCP receive buffer (default: `0`).

`aws.client.socketTimeout`
: The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out (default: `50000`).

`aws.client.storageClass`
: The S3 storage class applied to stored objects, one of \[`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`\] (default: `STANDARD`).

`aws.client.storageEncryption`
: The S3 server side encryption to be used when saving objects on S3. Can be `AES256` or `aws:kms` (default: none).

`aws.client.storageKmsKeyId`
: :::{versionadded} 22.05.0-edge
  :::
: The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket.

`aws.client.targetThroughputInGbps`
: :::{versionadded} 25.06.0-edge
  :::
: The target network throughput (in Gbps) used for S3 uploads and downloads (default: `10`).

`aws.client.userAgent`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The HTTP user agent header passed with all HTTP requests.

`aws.client.uploadChunkSize`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The size of a single part in a multipart upload (default: `100 MB`).

`aws.client.uploadMaxAttempts`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).

`aws.client.uploadMaxThreads`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The maximum number of threads used for multipart upload (default: `10`).

`aws.client.uploadRetrySleep`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The time to wait after a failed upload attempt to retry the part upload (default: `500ms`).

`aws.client.uploadStorageClass`
: The S3 storage class applied to stored objects. Can be `STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, or `INTELLIGENT_TIERING` (default: `STANDARD`).

(config-azure)=

## `azure`

The `azure` scope allows you to configure the interactions with Azure, including Azure Batch and Azure Blob Storage.

The following settings are available:

`azure.activeDirectory.servicePrincipalId`
: The service principal client ID. Defaults to environment variable `AZURE_CLIENT_ID`.

`azure.activeDirectory.servicePrincipalSecret`
: The service principal client secret. Defaults to environment variable `AZURE_CLIENT_SECRET`.

`azure.activeDirectory.tenantId`
: The Azure tenant ID. Defaults to environment variable `AZURE_TENANT_ID`.

`azure.azcopy.blobTier`
: The blob [access tier](https://learn.microsoft.com/en-us/azure/storage/blobs/access-tiers-overview) used by `azcopy` to upload files to Azure Blob Storage. Valid options are `None`, `Hot`, or `Cool` (default: `None`).

`azure.azcopy.blockSize`
: The block size (in MB) used by `azcopy` to transfer files between Azure Blob Storage and compute nodes (default: `4`).

`azure.batch.accountKey`
: The batch service account key. Defaults to environment variable `AZURE_BATCH_ACCOUNT_KEY`.

`azure.batch.accountName`
: The batch service account name. Defaults to environment variable `AZURE_BATCH_ACCOUNT_NAME`.

`azure.batch.allowPoolCreation`
: Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: `false`).

`azure.batch.autoPoolMode`
: Enable the automatic creation of batch pools depending on the pipeline resources demand (default: `true`).

`azure.batch.copyToolInstallMode`
: The mode in which the `azcopy` tool is installed by Nextflow (default: `'node'`).

: The following options are available:

  `'node'`
  : The `azcopy` tool is installed once during the pool creation.

  `'task'`
  : The `azcopy` tool is installed for each task execution.

  `'off'`
  : The `azcopy` tool is not installed.

`azure.batch.deleteJobsOnCompletion`
: :::{versionchanged} 23.08.0-edge
  Default value was changed from `true` to `false`.
  :::
: Delete all jobs when the workflow completes (default: `false`).

`azure.batch.deletePoolsOnCompletion`
: Delete all compute node pools when the workflow completes (default: `false`).

`azure.batch.deleteTasksOnCompletion`
: :::{versionadded} 23.08.0-edge
  :::
: Delete each task when it completes (default: `true`).
: Although this setting is enabled by default, failed tasks will not be deleted unless it is explicitly enabled. This way, the default behavior is that successful tasks are deleted while failed tasks are preserved for debugging purposes.

`azure.batch.endpoint`
: The batch service endpoint e.g. `https://nfbatch1.westeurope.batch.azure.com`.

`azure.batch.jobMaxWallClockTime`
: :::{versionadded} 25.04.0
  :::
: The maximum elapsed time that jobs may run, measured from the time they are created (default: `30d`).
: If jobs do not complete within this time limit, the Batch service terminates them and any tasks still running.

`azure.batch.location`
: The name of the batch service region, e.g. `westeurope` or `eastus2`. Not needed when the endpoint is specified.

`azure.batch.poolIdentityClientId`
: :::{versionadded} 25.05.0-edge
  :::
: The client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) that is available on all Azure Batch node pools. This identity is used by Fusion to authenticate to Azure storage. If set to `'auto'`, Fusion will use the first available managed identity.

`azure.batch.pools.<name>.autoScale`
: Enable autoscaling feature for the pool identified with `<name>`.

`azure.batch.pools.<name>.fileShareRootPath`
: The internal root mount point when mounting File Shares. Must be `/mnt/resource/batch/tasks/fsmounts` for CentOS nodes or `/mnt/batch/tasks/fsmounts` for Ubuntu nodes (default: CentOS).

`azure.batch.pools.<name>.lowPriority`
: Enable the use of low-priority VMs (default: `false`).

: :::{warning}
  As of September 30, 2025, Low Priority VMs will no longer be supported in Azure Batch accounts that use Batch Managed mode for pool allocation. You may continue to use this setting to configure Spot VMs in Batch accounts configured with User Subscription mode.
  :::

`azure.batch.pools.<name>.maxVmCount`
: The max number of virtual machines when using auto scaling.

`azure.batch.pools.<name>.mountOptions`
: The mount options for mounting the file shares (default: `-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp`).

`azure.batch.pools.<name>.offer`
: The offer type of the virtual machine type used by the pool identified with `<name>` (default: `centos-container`).

`azure.batch.pools.<name>.privileged`
: Enable the task to run with elevated access. Ignored if `runAs` is set (default: `false`).

`azure.batch.pools.<name>.publisher`
: The publisher of virtual machine type used by the pool identified with `<name>` (default: `microsoft-azure-batch`).

`azure.batch.pools.<name>.runAs`
: The username under which the task is run. The user must already exist on each node of the pool.

`azure.batch.pools.<name>.scaleFormula`
: The [scale formula](https://docs.microsoft.com/en-us/azure/batch/batch-automatic-scaling) for the pool identified with `<name>`.

`azure.batch.pools.<name>.scaleInterval`
: The interval at which to automatically adjust the Pool size according to the autoscale formula. Must be at least 5 minutes and at most 168 hours (default: `10 mins`).

`azure.batch.pools.<name>.schedulePolicy`
: The scheduling policy for the pool identified with `<name>`. Can be either `spread` or `pack` (default: `spread`).

`azure.batch.pools.<name>.sku`
: The ID of the Compute Node agent SKU which the pool identified with `<name>` supports (default: `batch.node.centos 8`).

`azure.batch.pools.<name>.startTask.privileged`
: :::{versionadded} 24.03.0-edge
  :::
: Enable the `startTask` to run with elevated access (default: `false`).

`azure.batch.pools.<name>.startTask.script`
: :::{versionadded} 24.03.0-edge
  :::
: The `startTask` that is executed as the node joins the Azure Batch node pool.

`azure.batch.pools.<name>.virtualNetwork`
: :::{versionadded} 23.03.0-edge
  :::
: The subnet ID of a virtual network in which to create the pool.

`azure.batch.pools.<name>.vmCount`
: The number of virtual machines provisioned by the pool identified with `<name>`.

`azure.batch.pools.<name>.vmType`
: The virtual machine type used by the pool identified with `<name>`.

`azure.batch.terminateJobsOnCompletion`
: :::{versionadded} 23.05.0-edge
  :::
: When the workflow completes, set all jobs to terminate on task completion (default: `true`).

`azure.managedIdentity.clientId`
: The client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview). Defaults to environment variable `AZURE_MANAGED_IDENTITY_USER`.

`azure.managedIdentity.system`
: When `true`, use the system-assigned [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) to authenticate Azure resources. Defaults to environment variable `AZURE_MANAGED_IDENTITY_SYSTEM`.

`azure.registry.password`
: The password to connect to a private container registry.

`azure.registry.server`
: The container registry from which to pull the Docker images (default: `docker.io`).

`azure.registry.userName`
: The username to connect to a private container registry.

`azure.retryPolicy.delay`
: Delay when retrying failed API requests (default: `250ms`).

`azure.retryPolicy.jitter`
: Jitter value when retrying failed API requests (default: `0.25`).

`azure.retryPolicy.maxAttempts`
: Max attempts when retrying failed API requests (default: `10`).

`azure.retryPolicy.maxDelay`
: Max delay when retrying failed API requests (default: `90s`).

`azure.storage.accountKey`
: The blob storage account key. Defaults to environment variable `AZURE_STORAGE_ACCOUNT_KEY`.

`azure.storage.accountName`
: The blob storage account name. Defaults to environment variable `AZURE_STORAGE_ACCOUNT_NAME`.

`azure.storage.fileShares.<name>.mountOptions`
: The file share mount options.

`azure.storage.fileShares.<name>.mountPath`
: The file share mount path.

`azure.storage.sasToken`
: The blob storage shared access signature (SAS) token, which can be provided instead of an account key. Defaults to environment variable `AZURE_STORAGE_SAS_TOKEN`.

`azure.storage.tokenDuration`
: The duration of the SAS token generated by Nextflow when the `sasToken` option is *not* specified (default: `48h`).

(config-charliecloud)=

## `charliecloud`

The `charliecloud` scope controls how [Charliecloud](https://hpc.github.io/charliecloud/) containers are executed by Nextflow.

The following settings are available:

`charliecloud.cacheDir`
: The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`charliecloud.enabled`
: Execute tasks with Charliecloud containers (default: `false`).

`charliecloud.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`charliecloud.pullTimeout`
: The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: `20 min`).

`charliecloud.registry`
: The registry from where images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`charliecloud.runOptions`
: Specify extra command line options supported by the `ch-run` command.

`charliecloud.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.

`charliecloud.writableInputMounts`
: When `false`, mount input directories as read-only (default: `true`).

`charliecloud.writeFake`
: Run containers from storage in writeable mode using overlayfs (default: `true`).
: This option requires unprivileged `overlayfs` (Linux kernel >= 5.11). For full support, tempfs with xattrs in the user namespace (Linux kernel >= 6.6) is required. See [charliecloud documentation](https://hpc.github.io/charliecloud/ch-run.html#ch-run-overlay) for details.

(config-conda)=

## `conda`

The `conda` scope controls the creation of Conda environments by the Conda package manager.

The following settings are available:

`conda.cacheDir`
: The path where Conda environments are stored. It should be accessible from all compute nodes when using a shared file system.

`conda.channels`
: The list of Conda channels that can be used to resolve Conda packages. Channel priority decreases from left to right.

`conda.createOptions`
: Extra command line options for the `conda create` command. See the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/commands/create.html) for more information.

`conda.createTimeout`
: The amount of time to wait for the Conda environment to be created before failing (default: `20 min`).

`conda.enabled`
: Execute tasks with Conda environments (default: `false`).

`conda.useMamba`
: Use [Mamba](https://github.com/mamba-org/mamba) instead of `conda` to create Conda environments (default: `false`).

`conda.useMicromamba`
: :::{versionadded} 22.05.0-edge
  :::
: Use [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) instead of `conda` to create Conda environments (default: `false`).

(config-dag)=

## `dag`

The `dag` scope controls the workflow diagram generated by Nextflow.

The following settings are available:

`dag.depth`
: :::{versionadded} 23.10.0
  :::
: *Supported by the HTML and Mermaid renderers.*
: Controls the maximum depth at which to render sub-workflows (default: no limit).

`dag.direction`
: :::{versionadded} 23.10.0
  :::
: *Supported by the Graphviz, DOT, HTML and Mermaid renderers.*
: Controls the direction of the DAG, can be `'LR'` (left-to-right) or `'TB'` (top-to-bottom) (default: `'TB'`).

`dag.enabled`
: When `true` enables the generation of the DAG file (default: `false`).

`dag.file`
: Graph file name (default: `'dag-<timestamp>.html'`).

`dag.overwrite`
: When `true` overwrites any existing DAG file with the same name (default: `false`).

`dag.verbose`
: :::{versionadded} 23.10.0
  :::
: *Only supported by the HTML and Mermaid renderers.*
: When `false`, channel names are omitted, operators are collapsed, and empty workflow inputs are removed (default: `false`).

(config-docker)=

## `docker`

The `docker` scope controls how [Docker](https://www.docker.com) containers are executed by Nextflow.

The following settings are available:

`docker.enabled`
: Enable Docker execution (default: `false`).

`docker.engineOptions`
: Specify additional options supported by the Docker engine i.e. `docker [OPTIONS]`.

`docker.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`docker.fixOwnership`
: Fix ownership of files created by the Docker container (default: `false`).

`docker.legacy`
: Use command line options removed since Docker 1.10.0 (default: `false`).

`docker.mountFlags`
: Add the specified flags to the volume mounts e.g. `'ro,Z'`.

`docker.registry`
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`docker.registryOverride`
: :::{versionadded} 25.06.0-edge
  :::
: When `true`, forces the override of the registry name in fully qualified container image names with the registry specified by `docker.registry` (default: `false`).
: This setting allows you to redirect container image pulls from their original registry to a different registry, such as a private mirror or proxy.

`docker.remove`
: Clean up the container after the execution (default: `true`). See the [Docker documentation](https://docs.docker.com/engine/reference/run/#clean-up---rm) for details.

`docker.runOptions`
: Specify extra command line options supported by the `docker run` command. See the [Docker documentation](https://docs.docker.com/engine/reference/run/) for details.

`docker.sudo`
: Executes Docker run command as `sudo` (default: `false`).

`docker.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.

`docker.tty`
: Allocates a pseudo-tty (default: `false`).

`docker.writableInputMounts`
: When `false`, mount input directories as read-only (default: `true`).

(config-env)=

## `env`

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

## `executor`

The `executor` scope controls various executor behaviors.

The following settings are available:

`executor.account`
: :::{versionadded} 24.04.0
  :::
: *Used only by the {ref}`slurm-executor`, {ref}`lsf-executor`, {ref}`pbs-executor` and {ref}`pbspro-executor` executors.*
: The project or organization account that should be charged for running the pipeline jobs.

`executor.cpus`
: *Used only by the local executor.*
: The maximum number of CPUs made available by the underlying system.

`executor.dumpInterval`
: Determines how often to log the executor status (default: `5 min`).

`executor.exitReadTimeout`
: *Used only by grid executors.*
: Determines how long to wait for the `.exitcode` file to be created after the task has completed, before returning an error status (default: `270 sec`).

`executor.jobName`
: *Used only by grid executors and Google Batch.*
: Determines the name of jobs submitted to the underlying cluster executor:
  ```groovy
  executor.jobName = { "$task.name - $task.hash" }
  ```
: The job name should satisfy the validation constraints of the underlying scheduler.

`executor.killBatchSize`
: Determines the number of jobs that can be killed in a single command execution (default: `100`).

`executor.memory`
: *Used only by the local executor.*
: The maximum amount of memory made available by the underlying system.

`executor.name`
: The name of the executor to be used (default: `local`).

`executor.perCpuMemAllocation`
: :::{versionadded} 23.07.0-edge
  :::
: *Used only by the {ref}`slurm-executor` executor.*
: When `true`, memory allocations for SLURM jobs are specified as `--mem-per-cpu <task.memory / task.cpus>` instead of `--mem <task.memory>`.

`executor.perJobMemLimit`
: *Used only by the {ref}`lsf-executor` executor.*
: Enables the *per-job* memory limit mode for LSF jobs.

`executor.perTaskReserve`
: *Used only by the {ref}`lsf-executor` executor.*
: Enables the *per-task* memory reserve mode for LSF jobs.

`executor.pollInterval`
: Determines how often to check for process termination. Default varies for each executor.

`executor.queueGlobalStatus`
: :::{versionadded} 23.01.0-edge
  :::
: Determines how job status is retrieved. When `false` only the queue associated with the job execution is queried. When `true` the job status is queried globally i.e. irrespective of the submission queue (default: `false`).

`executor.queueSize`
: The number of tasks the executor will handle in a parallel manner. A queue size of zero corresponds to no limit. Default varies for each executor.

`executor.queueStatInterval`
: *Used only by grid executors.*
: Determines how often to fetch the queue status from the scheduler (default: `1 min`).

`executor.submitRateLimit`
: Determines the max rate of job submission per time unit, for example `'10sec'` (10 jobs per second) or `'50/2min'` (50 jobs every 2 minutes) (default: unlimited).

`executor.retry.delay`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Delay when retrying failed job submissions (default: `500ms`).

`executor.retry.jitter`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Jitter value when retrying failed job submissions (default: `0.25`).

`executor.retry.maxAttempts`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Max attempts when retrying failed job submissions (default: `3`).

`executor.retry.maxDelay`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Max delay when retrying failed job submissions (default: `30s`).

`executor.retry.reason`
: :::{versionadded} 22.03.0-edge
  :::
: :::{versionchanged} 25.10.0
  This option was renamed from `executor.submit.retry.reason` to `executor.retry.reason`.
  :::
: *Used only by grid executors.*
: Regex pattern that when verified causes a failed submit operation to be re-tried (default: `Socket timed out`).

### Executor-specific defaults

Some executor settings have different default values depending on the executor.

| Executor       | `queueSize` | `pollInterval` |
| -------------- | ----------- | -------------- |
| AWS Batch      | `1000`      | `10s`          |
| Azure Batch    | `1000`      | `10s`          |
| Google Batch   | `1000`      | `10s`          |
| Grid Executors | `100`       | `5s`           |
| Kubernetes     | `100`       | `5s`           |
| Local          | N/A         | `100ms`        |

### Executor-specific configuration

Executor config settings can be applied to specific executors by prefixing the executor name with the symbol `$` and using it as special scope. For example:

```groovy
// block syntax
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

// dot syntax
executor.$sge.queueSize = 100
executor.$sge.pollInterval = '30sec'
executor.$local.cpus = 8
executor.$local.memory = '32 GB'
```

(config-fusion)=

## `fusion`

The `fusion` scope provides advanced configuration for the use of the {ref}`Fusion file system <fusion-page>`.

The following settings are available:

`fusion.cacheSize`
: :::{versionadded} 23.11.0-edge
  :::
: The maximum size of the local cache used by the Fusion client.

`fusion.containerConfigUrl`
: The URL of the container layer that provides the Fusion client.

`fusion.enabled`
: Enable the Fusion file system (default: `false`).

`fusion.exportStorageCredentials`
: Export access credentials required by the underlying object storage as environment variables (e.g., `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, and `AWS_SESSION_TOKEN` for AWS S3) to task execution environments (default: `false`).

  :::{note}
  This configuration does not mount or provide access to credential files. For example, AWS credentials like `~/.aws/credentials`, `~/.aws/config`, and SSO cache files are not mounted. AWS SSO users must export credentials to environment variables:
  
  ```bash
  eval "$(aws configure export-credentials --format env)"
  ```
  :::

  :::{warning}
  This option leaks credentials is the task launcher script. It should only be used for testing and development purposes.
  :::

`fusion.logLevel`
: The log level of the Fusion client.

`fusion.logOutput`
: The output location of the Fusion log.

`fusion.privileged`
: :::{versionadded} 23.10.0
  :::
: Enable privileged containers for Fusion (default: `true`).
: Non-privileged use is supported only on Kubernetes with the [k8s-fuse-plugin](https://github.com/nextflow-io/k8s-fuse-plugin) or a similar FUSE device plugin.

`fusion.snapshots`
: :::{versionadded} 25.03.0-edge
  :::
: *Currently only supported for AWS Batch*
: Enable Fusion snapshotting (preview, default: `false`). This feature allows Fusion to automatically restore a job when it is interrupted by a spot reclamation.

`fusion.tags`
: *Currently only supported for S3.*
: The pattern that determines how tags are applied to files created via the Fusion client (default: `[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)`). Set to `false` to disable tags.

(config-google)=

## `google`

The `google` scope allows you to configure the interactions with Google Cloud, including Google Cloud Batch and Google Cloud Storage.

The following settings are available:

`google.enableRequesterPaysBuckets`
: Use the given Google Cloud project ID as the billing project for storage access (default: `false`). Required when accessing data from [requester pays](https://cloud.google.com/storage/docs/requester-pays) buckets.

`google.httpConnectTimeout`
: :::{versionadded} 23.06.0-edge
  :::
: The HTTP connection timeout for Cloud Storage API requests (default: `'60s'`).

`google.httpReadTimeout`
: :::{versionadded} 23.06.0-edge
  :::
: The HTTP read timeout for Cloud Storage API requests (default: `'60s'`).

`google.location`
: The Google Cloud location where jobs are executed (default: `us-central1`).

`google.project`
: The Google Cloud project ID to use for pipeline execution.

`google.batch.allowedLocations`
: :::{versionadded} 22.12.0-edge
  :::
: The set of [allowed locations](https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy) for VMs to be provisioned (default: no restriction).

`google.batch.autoRetryExitCodes`
: :::{versionadded} 24.07.0-edge
  :::
: The list of exit codes that should be automatically retried by Google Batch when `google.batch.maxSpotAttempts` is greater than 0 (default: `[50001]`).
: See [Google Batch documentation](https://cloud.google.com/batch/docs/troubleshooting#reserved-exit-codes) for the complete list of retryable exit codes.

`google.batch.bootDiskImage`
: :::{versionadded} 24.08.0-edge
  :::
: The image URI of the virtual machine boot disk, e.g `batch-debian` (default: none).
: See [Google documentation](https://cloud.google.com/batch/docs/vm-os-environment-overview#vm-os-image-options) for details.

`google.batch.bootDiskSize`
: The size of the virtual machine boot disk, e.g `50.GB` (default: none).

`google.batch.cpuPlatform`
: The [minimum CPU Platform](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications), e.g. `'Intel Skylake'` (default: none).

`google.batch.gcsfuseOptions`
: :::{versionadded} 25.03.0-edge
  :::
: List of custom mount options for `gcsfuse` (default: `['-o rw', '-implicit-dirs']`).

`google.batch.installOpsAgent`
: Enable the installation of the Ops Agent on Google Batch instances for enhanced monitoring and logging (default: `false`).

: :::{note}
  The Ops Agent requires a compatible boot disk image. For Google Batch, use [Batch-debian images](https://docs.cloud.google.com/batch/docs/vm-os-environment-overview#vm-os-image-options) (e.g., `batch-debian`) with `google.batch.bootDiskImage`. The default Container-Optimized OS (`batch-cos`) is not compatible with the Ops Agent.
  :::

`google.batch.logsPath`
: :::{versionadded} 25.11.0-edge
  :::
: The Google Cloud Storage path where job logs should be stored, e.g. `gs://my-logs-bucket/logs`.
: When specified, Google Batch will write job logs to this location instead of [Cloud Logging](https://cloud.google.com/logging/docs). The bucket must be accessible and writable by the service account.

`google.batch.maxSpotAttempts`
: :::{versionadded} 23.11.0-edge
  :::
: :::{versionchanged} 24.08.0-edge
  The default value was changed from `5` to `0`.
  :::
: Max number of execution attempts of a job interrupted by a Compute Engine Spot reclaim event (default: `0`).
: See also: `google.batch.autoRetryExitCodes`

`google.batch.network`
: The URL of an existing network resource to which the VM will be attached.

: You can specify the network as a full or partial URL. For example, the following are all valid URLs:

  - `https://www.googleapis.com/compute/v1/projects/{project}/global/networks/{network}`
  - `projects/{project}/global/networks/{network}`
  - `global/networks/{network}`

`google.batch.networkTags`
: The [network tags](https://cloud.google.com/vpc/docs/add-remove-network-tags) to be applied to the instances created by Google Batch jobs (e.g., `['allow-ssh', 'allow-http']`).
: Network tags are ignored when using instance templates.

`google.batch.serviceAccountEmail`
: The Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.
: This service account will only be used for tasks submitted by Nextflow, not for Nextflow itself. See {ref}`google-credentials` for more information on Google Cloud credentials.

`google.batch.spot`
: Enable the use of spot virtual machines (default: `false`).

`google.batch.subnetwork`
: The URL of an existing subnetwork resource in the network to which the VM will be attached.

: You can specify the subnetwork as a full or partial URL. For example, the following are all valid URLs:

  - `https://www.googleapis.com/compute/v1/projects/{project}/regions/{region}/subnetworks/{subnetwork}`
  - `projects/{project}/regions/{region}/subnetworks/{subnetwork}`
  - `regions/{region}/subnetworks/{subnetwork}`

`google.batch.usePrivateAddress`
: Do not provision public IP addresses for VMs, such that they only have an internal IP address (default: `false`).
: When this option is enabled, jobs can only load Docker images from Google Container Registry, and cannot use external services other than Google APIs.

`google.storage.retryPolicy.maxAttempts`
: :::{versionadded} 23.11.0-edge
  :::
: Max attempts when retrying failed API requests to Cloud Storage (default: `10`).

`google.storage.retryPolicy.maxDelay`
: :::{versionadded} 23.11.0-edge
  :::
: Max delay when retrying failed API requests to Cloud Storage (default: `'90s'`).

`google.storage.retryPolicy.multiplier`
: :::{versionadded} 23.11.0-edge
  :::
: Delay multiplier when retrying failed API requests to Cloud Storage (default: `2.0`).

(config-k8s)=

## `k8s`

The `k8s` scope controls the deployment and execution of workflow applications in a Kubernetes cluster.

The following settings are available:

`k8s.autoMountHostPaths`
: Automatically mount host paths into the task pods (default: `false`). Only intended for development purposes when using a single node.

`k8s.cleanup`
: When `true`, successful pods are automatically deleted (default: `true`).

`k8s.client`
: Map of options for the Kubernetes HTTP client.
: If this option is specified, it will be used instead of `.kube/config`.
: The following options are available:
  - `server`
  - `token`
  - `tokenFile`
  - `verifySsl`
  - `sslCert`
  - `sslCertFile`
  - `clientCert`
  - `clientCertFile`
  - `clientKey`
  - `clientKeyFile`

`k8s.computeResourceType`
: :::{versionadded} 22.05.0-edge
  :::
: Whether to use Kubernetes `Pod` or `Job` resource type to carry out Nextflow tasks (default: `Pod`).

`k8s.context`
: The Kubernetes [configuration context](https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/) to use.

`k8s.cpuLimits`
: :::{versionadded} 24.04.0
  :::
: When `true`, set both the pod CPUs `request` and `limit` to the value specified by the `cpus` directive, otherwise set only the `request` (default: `false`).
: This setting is useful when a K8s cluster requires a CPU limit to be defined through a [LimitRange](https://kubernetes.io/docs/concepts/policy/limit-range/).

`k8s.fetchNodeName`
: :::{versionadded} 22.05.0-edge
  :::
: Include the hostname of each task in the execution trace (default: `false`).

`k8s.fuseDevicePlugin`
: :::{versionadded} 24.01.0-edge
  :::
: The FUSE device plugin to be used when enabling Fusion in unprivileged mode (default: `['nextflow.io/fuse': 1]`).

`k8s.httpConnectTimeout`
: :::{versionadded} 22.10.0
  :::
: The Kubernetes HTTP client request connection timeout e.g. `'60s'`.

`k8s.httpReadTimeout`
: :::{versionadded} 22.10.0
  :::
: The Kubernetes HTTP client request connection read timeout e.g. `'60s'`.

`k8s.imagePullPolicy`
: The strategy for pulling container images. Can be `IfNotPresent`, `Always`, `Never`.

`k8s.launchDir`
: The path where the workflow is launched and the user data is stored (default: `<volume-claim-mount-path>/<user-name>`). Must be a path in a shared K8s persistent volume.

`k8s.namespace`
: The Kubernetes namespace to use (default: `default`).

`k8s.pod`
: Additional pod configuration options such as environment variables, config maps, secrets, etc. Allows the same settings as the {ref}`process-pod` process directive.
: When using the `kuberun` command, this setting also applies to the submitter pod.

`k8s.projectDir`
: The path where Nextflow projects are downloaded (default: `<volume-claim-mount-path>/projects`). Must be a path in a shared K8s persistent volume.

`k8s.runAsUser`
: The user ID to be used to run the containers. Shortcut for the `securityContext` option.

`k8s.securityContext`
: The [security context](https://kubernetes.io/docs/tasks/configure-pod-container/security-context/) to use for all pods.

`k8s.serviceAccount`
: The Kubernetes [service account name](https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/) to use.

`k8s.storageClaimName`
: The name of the persistent volume claim where the shared work directory is stored.

`k8s.storageMountPath`
: The mount path for the persistent volume claim (default: `/workspace`).

`k8s.storageSubPath`
: The path in the persistent volume to be mounted (default: `/`).

`k8s.workDir`
: The path of the shared work directory (default: `<user-dir>/work`). Must be a path in a shared K8s persistent volume.

`k8s.debug.yaml`
: Save the pod spec for each task to `.command.yaml` in the task directory (default: `false`).

`k8s.retryPolicy.delay`
: Delay when retrying failed API requests (default: `250ms`).

`k8s.retryPolicy.jitter`
: Jitter value when retrying failed API requests (default: `0.25`).

`k8s.retryPolicy.maxAttempts`
: Max attempts when retrying failed API requests (default: `4`).

`k8s.retryPolicy.maxDelay`
: Max delay when retrying failed API requests (default: `90s`).

(config-lineage)=

## `lineage`

The `lineage` scope controls the generation of {ref}`cli-lineage` metadata.

The following settings are available:

`lineage.enabled`
: Enable generation of lineage metadata (default: `false`).

`lineage.store.location`
: The location of the lineage metadata store (default: `./.lineage`).

(config-mail)=

## `mail`

The `mail` scope controls the mail server used to send email notifications.

The following settings are available:

`mail.debug`
: Enable Java Mail logging for debugging purposes (default: `false`).

`mail.from`
: Default email sender address.

`mail.smtp.host`
: Host name of the mail server.

`mail.smtp.password`
: User password to connect to the mail server.

`mail.smtp.port`
: Port number of the mail server.

`mail.smtp.user`
: User name to connect to the mail server.

`mail.smtp.proxy.host`
: Host name of an HTTP web proxy server that will be used for connections to the mail server.

`mail.smtp.proxy.port`
: Port number for the HTTP web proxy server.

`mail.smtp.*`
: Any SMTP configuration property supported by the [Java Mail API](https://javaee.github.io/javamail/), which Nextflow uses to send emails. See the table of available properties [here](https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties).

(config-manifest)=

## `manifest`

The `manifest` scope allows you to define some metadata that is useful when publishing or running your pipeline.

The following settings are available:

`manifest.author`
: :::{deprecated} 24.09.0-edge
  Use `manifest.contributors` instead.
  :::
: Project author name (use a comma to separate multiple names).

`manifest.contributors`
: :::{versionadded} 24.09.0-edge
  :::
: List of project contributors. Should be a list of maps.

: The following fields are supported in the contributor map:

  `name`
  : The contributor name.

  `affiliation`
  : The contributor affiliated organization.

  `email`
  : The contributor email address.

  `github`
  : The contributor GitHub URL.

  `contribution`
  : List of contribution types, each element can be one of `'author'`, `'maintainer'`, or `'contributor'`.

  `orcid`
  : The contributor [ORCID](https://orcid.org/) URL.

`manifest.defaultBranch`
: Git repository default branch (default: `master`).

`manifest.description`
: Free text describing the workflow project.

`manifest.docsUrl`
: Project documentation URL.

`manifest.doi`
: Project related publication DOI identifier.

`manifest.gitmodules`
: Controls whether git sub-modules should be cloned with the main repository.
: Can be either a boolean value, a list of submodule names, or a comma-separated string of submodule names.

`manifest.homePage`
: Project home page URL.

`manifest.icon`
: Project related icon location (Relative path or URL).

`manifest.license`
: Project license.

`manifest.mainScript`
: Project main script (default: `main.nf`).

`manifest.name`
: Project short name.

`manifest.nextflowVersion`
: Minimum required Nextflow version.

: This setting may be useful to ensure that a specific version is used:

  ```groovy
  manifest.nextflowVersion = '1.2.3'        // exact match
  manifest.nextflowVersion = '1.2+'         // 1.2 or later (excluding 2 and later)
  manifest.nextflowVersion = '>=1.2'        // 1.2 or later
  manifest.nextflowVersion = '>=1.2, <=1.5' // any version in the 1.2 .. 1.5 range
  manifest.nextflowVersion = '!>=1.2'       // with ! prefix, stop execution if current version does not match required version.
  ```

: See {ref}`stdlib-types-versionnumber` for details.

`manifest.organization`
: Project organization.

`manifest.recurseSubmodules`
: Pull submodules recursively from the Git repository.

`manifest.version`
: Project version number.

(config-nextflow)=

## `nextflow`

:::{versionchanged} 24.10.0
The `nextflow.publish.retryPolicy` settings were moved to `workflow.output.retryPolicy`.
:::

:::{versionchanged} 25.06.0-edge
The `workflow.output.retryPolicy` settings were moved to `nextflow.retryPolicy`.
:::

`retryPolicy.delay`
: Delay used for retryable operations (default: `350ms`).

`retryPolicy.jitter`
: Jitter value used for retryable operations (default: `0.25`).

`retryPolicy.maxAttempts`
: Max attempts used for retryable operations (default: `5`).

`retryPolicy.maxDelay`
: Max delay used for retryable operations (default: `90s`).

(config-notification)=

## `notification`

The `notification` scope controls the automatic sending of an email notification on workflow completion.

The following settings are available:

`notification.attributes`
: Map of variables that can be used in the template file.

`notification.enabled`
: Send an email notification when the workflow execution completes (default: `false`).

`notification.from`
: Sender address for the email notification.

`notification.template`
: Path of a template file containing the contents of the email notification.

`notification.to`
: Recipient address for the email notification. Multiple addresses can be specified as a comma-separated list.

(config-podman)=

## `podman`

The `podman` scope controls how [Podman](https://podman.io/) containers are executed by Nextflow.

The following settings are available:

`podman.enabled`
: Execute tasks with Podman containers (default: `false`).

`podman.engineOptions`
: Specify additional options supported by the Podman engine i.e. `podman [OPTIONS]`.

`podman.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`podman.mountFlags`
: Add the specified flags to the volume mounts e.g. `'ro,Z'`.

`podman.registry`
: The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`podman.remove`
: Clean-up the container after the execution (default: `true`).

`podman.runOptions`
: Specify extra command line options supported by the `podman run` command.

`podman.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `'auto'` to create a temporary directory each time a container is created.

(config-report)=

## `report`

The `report` scope allows you to configure the workflow {ref}`execution-report`.

The following settings are available:

`report.enabled`
: Create the execution report on workflow completion (default: `false`).

`report.file`
: The path of the created execution report file (default: `'report-<timestamp>.html'`).

`report.overwrite`
: Overwrite any existing report file with the same name (default: `false`).

(config-sarus)=

## `sarus`

The `sarus` scope controls how [Sarus](https://sarus.readthedocs.io) containers are executed by Nextflow.

The following settings are available:

`sarus.enabled`
: Execute tasks with Sarus containers (default: `false`).

`sarus.envWhitelist`
: Comma-separated list of environment variable names to be included in the container environment.

`sarus.runOptions`
: Specify extra command line options supported by the `sarus run` command.
: See the [Sarus user guide](https://sarus.readthedocs.io/en/stable/user/user_guide.html) for details.

`sarus.tty`
: Allocates a pseudo-tty (default: `false`).

(config-shifter)=

## `shifter`

The `shifter` scope controls how [Shifter](https://docs.nersc.gov/programming/shifter/overview/) containers are executed by Nextflow.

The following settings are available:

`shifter.enabled`
: Execute tasks with Shifter containers (default: `false`).

`shifter.envWhitelist`
: Comma-separated list of environment variable names to be included in the container environment.

(config-singularity)=

## `singularity`

The `singularity` scope controls how [Singularity](https://sylabs.io/singularity/) containers are executed by Nextflow.

The following settings are available:

`singularity.autoMounts`
: :::{versionchanged} 23.09.0-edge
  Default value was changed from `false` to `true`.
  :::
: Automatically mount host paths in the executed container (default: `true`). It requires the `user bind control` feature to be enabled in your Singularity installation.

`singularity.cacheDir`
: The directory where remote Singularity images are stored. When using a compute cluster, it must be a shared folder accessible to all compute nodes.

`singularity.enabled`
: Execute tasks with Singularity containers (default: `false`).

`singularity.engineOptions`
: Specify additional options supported by the Singularity engine i.e. `singularity [OPTIONS]`.

`singularity.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`singularity.libraryDir`
: Directory where remote Singularity images are retrieved. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`singularity.noHttps`
: Pull the Singularity image with http protocol (default: `false`).

`singularity.ociAutoPull`
: :::{versionadded} 23.12.0-edge
  :::
: *Requires Singularity 3.11 or later*
: When enabled, OCI (and Docker) container images are pull and converted to a SIF image file format implicitly by the Singularity run command, instead of Nextflow (default: `false`).

  :::{note}
  Leave `ociAutoPull` disabled if willing to build a Singularity native image with Wave (see the {ref}`wave-singularity` section).
  :::

`singularity.ociMode`
: :::{versionadded} 23.12.0-edge
  :::
: *Requires Singularity 4 or later*
: Enable OCI-mode, that allows running native OCI compliant container image with Singularity using `crun` or `runc` as low-level runtime (default: `false`).
: See `--oci` flag in the [Singularity documentation](https://docs.sylabs.io/guides/4.0/user-guide/oci_runtime.html#oci-mode) for more details and requirements (default: `false`).

  :::{note}
  Leave `ociMode` disabled if you are willing to build a Singularity native image with Wave (see the {ref}`wave-singularity` section).
  :::

`singularity.pullTimeout`
: The amount of time the Singularity pull can last, after which the process is terminated (default: `20 min`).

`singularity.registry`
: :::{versionadded} 22.12.0-edge
  :::
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`singularity.runOptions`
: Specify extra command line options supported by `singularity exec`.

(config-spack)=

## `spack`

The `spack` scope controls the creation of a Spack environment by the Spack package manager.

The following settings are available:

`spack.cacheDir`
: The path where Spack environments are stored. It should be accessible from all compute nodes when using a shared file system.

`spack.checksum`
: Enable checksum verification of source tarballs (default: `true`).
: Only disable when requesting a package version not yet encoded in the corresponding Spack recipe.

`spack.createTimeout`
: The amount of time to wait for the Spack environment to be created before failing (default: `60 min`).

`spack.enabled`
: Execute tasks with Spack environments (default: `false`).

`spack.parallelBuilds`
: The maximum number of parallel package builds (default: the number of available CPUs).

(config-timeline)=

## `timeline`

The `timeline` scope controls the execution timeline report generated by Nextflow.

The following settings are available:

`timeline.enabled`
: Create the timeline report on workflow completion file (default: `false`).

`timeline.file`
: Timeline file name (default: `'timeline-<timestamp>.html'`).

`timeline.overwrite`
: Overwrite any existing timeline file with the same name (default: `false`).

(config-tower)=

## `tower`

The `tower` scope controls the settings for [Seqera Platform](https://seqera.io) (formerly Tower Cloud).

The following settings are available:

`tower.accessToken`
: The unique access token for your Seqera Platform account.
: Your `accessToken` can be obtained from your Seqera Platform instance in the [Tokens page](https://cloud.seqera.io/tokens).

`tower.computeEnvId`
: The compute environment ID in your Seqera Platform account used to launch pipelines (default: the primary compute environment in the selected workspace).

`tower.enabled`
: Enable workflow monitoring with Seqera Platform (default: `false`).

`tower.endpoint`
: The endpoint of your Seqera Platform instance (default: `https://api.cloud.seqera.io`).

`tower.workspaceId`
: The workspace ID in Seqera Platform in which to save the run (default: the launching user's personal workspace).
: The workspace ID can also be specified using the environment variable `TOWER_WORKSPACE_ID` (config file has priority over the environment variable).

(config-trace)=

## `trace`

The `trace` scope controls the layout of the execution trace file generated by Nextflow.

The following settings are available:

`trace.enabled`
: Create the execution trace file on workflow completion (default: `false`).

`trace.fields`
: Comma-separated list of {ref}`trace fields <trace-fields>` to include in the report.

`trace.file`
: Trace file name (default: `'trace-<timestamp>.txt'`).

`trace.overwrite`
: Overwrite any existing trace file with the same name (default: `false`).

`trace.raw`
: Report trace metrics as raw numbers where applicable, i.e. report duration values in milliseconds and memory values in bytes (default: `false`).

`trace.sep`
: Character used to separate values in each row (default: `\t`).

(config-wave)=

## `wave`

The `wave` scope provides advanced configuration for the use of {ref}`Wave containers <wave-page>`.

The following settings are available:

`wave.enabled`
: Enable the use of Wave containers (default: `false`).

`wave.endpoint`
: The Wave service endpoint (default: `https://wave.seqera.io`).

`wave.freeze`
: :::{versionadded} 23.07.0-edge
  :::
: Enable Wave container freezing (default: `false`). Wave will provision a non-ephemeral container image that will be pushed to a container repository of your choice.
: The target registry must be specified using the `wave.build.repository` setting. It is also recommended to specify a custom cache repository using `wave.build.cacheRepository`.
: :::{note}
  The container registry authentication must be managed by the underlying infrastructure.
  :::

`wave.mirror`
: :::{versionadded} 24.09.1-edge
  :::
: Enable Wave container mirroring (default: `false`). Wave will mirror (i.e. copy) the containers in your pipeline to a container registry of your choice, so that pipeline tasks can pull the containers from this registry instead of the original one.
: The mirrored containers will have the same name, digest, and metadata.
: The target registry must be specified using the `wave.build.repository` setting.
: This option is only compatible with `wave.strategy = 'container'`. It cannot be used with `wave.freeze`.
: :::{note}
  The container registry authentication must be managed by the underlying infrastructure.
  :::

`wave.strategy`
: The strategy to be used when resolving multiple Wave container requirements (default: `'container,dockerfile,conda'`).

`wave.build.cacheRepository`
: The container repository used to cache image layers built by the Wave service.
: The corresponding credentials must be provided in your Seqera Platform account.

`wave.build.compression.mode`
: :::{versionadded} 25.05.0-edge
  :::
: The compression algorithm that should be used when building the container. Allowed values are: `gzip`, `estargz` and `zstd` (default: `gzip`).

`wave.build.compression.level`
: :::{versionadded} 25.05.0-edge
  :::
: Level of compression used when building a container depending the chosen algorithm: gzip, estargz (0-9) and zstd (0-22).

`wave.build.compression.force`
: :::{versionadded} 25.05.0-edge
  :::
: Forcefully apply compression option to all layers, including already existing layers (default: `false`).

`wave.build.conda.basePackages`
: One or more Conda packages to be always added in the resulting container (default: `conda-forge::procps-ng`).

`wave.build.conda.commands`
: One or more commands to be added to the Dockerfile used to build a Conda based image.

`wave.build.conda.mambaImage`
: The Mamba container image is used to build Conda based container. This is expected to be [micromamba-docker](https://github.com/mamba-org/micromamba-docker) image.

`wave.build.repository`
: The container repository where images built by Wave are uploaded.
: The corresponding credentials must be provided in your Seqera Platform account.

`wave.httpClient.connectTimeout`
: :::{versionadded} 22.06.0-edge
  :::
: The connection timeout for the Wave HTTP client (default: `30s`).

`wave.httpClient.maxRate`
: :::{versionadded} 25.01.0-edge
  :::
: The maximum request rate for the Wave HTTP client (default: `1/sec`).

`wave.retryPolicy.delay`
: :::{versionadded} 22.06.0-edge
  :::
: The initial delay when a failing HTTP request is retried (default: `450ms`).

`wave.retryPolicy.jitter`
: :::{versionadded} 22.06.0-edge
  :::
: The jitter factor used to randomly vary retry delays (default: `0.25`).

`wave.retryPolicy.maxAttempts`
: :::{versionadded} 22.06.0-edge
  :::
: The max number of attempts a failing HTTP request is retried (default: `5`).

`wave.retryPolicy.maxDelay`
: :::{versionadded} 22.06.0-edge
  :::
: The max delay when a failing HTTP request is retried (default: `90s`).

`wave.scan.allowedLevels`
: :::{versionadded} 24.09.1-edge
  :::
: Comma-separated list of allowed vulnerability levels when scanning containers for security vulnerabilities in `required` mode.

: Allowed values are: `low`, `medium`, `high`, `critical`.

: This option requires `wave.scan.mode = 'required'`.

`wave.scan.mode`
: :::{versionadded} 24.09.1-edge
  :::
: Enable Wave container security scanning. Wave will scan the containers in your pipeline for security vulnerabilities.

: The following options can be specified:

  `'none'`
  : No security scanning is performed.

  `'async'`
  : The containers used by your pipeline are scanned for security vulnerabilities. The task execution is carried out regardless of the security scan result.

  `'required'`
  : The containers used by your pipeline are scanned for security vulnerabilities. The task is only executed if the corresponding container is free of vulnerabilities.

(config-workflow)=

## `workflow`

:::{versionadded} 24.10.0
:::

The `workflow` scope provides workflow execution options.

The following settings are available:

`workflow.failOnIgnore`
: :::{versionadded} 24.05.0-edge
  :::
: When `true`, the pipeline will exit with a non-zero exit code if any failed tasks are ignored using the `ignore` {ref}`error strategy <process-error-strategy>` (default: `false`).

`workflow.onComplete`
: Specify a closure that will be invoked at the end of a workflow run (including failed runs). See {ref}`workflow-handlers` for more information.

`workflow.onError`
: Specify a closure that will be invoked if a workflow run is terminated. See {ref}`workflow-handlers` for more information.

`workflow.output.contentType`
: *Currently only supported for S3.*
: Specify the media type, also known as [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/MIME_types), of published files (default: `false`). Can be a string (e.g. `'text/html'`), or `true` to infer the content type from the file extension.

`workflow.output.copyAttributes`
: :::{versionadded} 25.01.0-edge
  :::
: *Currently only supported for local and shared filesystems.*
: Copy file attributes (such as the last modified timestamp) to the published file (default: `false`).

`workflow.output.enabled`
: Enable or disable publishing (default: `true`).

`workflow.output.ignoreErrors`
: When `true`, the workflow will not fail if a file can't be published for some reason (default: `false`).

`workflow.output.mode`
: The file publishing method (default: `'symlink'`).

: The following options are available:

  `'copy'`
  : Copy each file into the output directory.

  `'copyNoFollow'`
  : Copy each file into the output directory without following symlinks, i.e. only the link is copied.

  `'link'`
  : Create a hard link in the output directory for each file.

  `'move'`
  : Move each file into the output directory.
  : Should only be used for files which are not used by downstream processes in the workflow.

  `'rellink'`
  : Create a relative symbolic link in the output directory for each file.

  `'symlink'`
  : Create an absolute symbolic link in the output directory for each output file.

`workflow.output.overwrite`
: When `true` any existing file in the specified folder will be overwritten (default: `'standard'`).

: The following options are available:

  `false`
  : Never overwrite existing files.

  `true`
  : Always overwrite existing files.

  `'deep'`
  : Overwrite existing files when the file content is different.

  `'lenient'`
  : Overwrite existing files when the file size is different.

  `'standard'`
  : Overwrite existing files when the file size or last modified timestamp is different.

`workflow.output.storageClass`
: *Currently only supported for S3.*
: Specify the storage class for published files.

`workflow.output.tags`
: *Currently only supported for S3.*
: Specify arbitrary tags for published files.

: For example:
  ```groovy
  workflow.output.tags = [FOO: 'hello', BAR: 'world']
  ```
