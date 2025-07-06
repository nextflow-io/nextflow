(config-options)=

# Configuration options

This page lists all of the available settings in the {ref}`Nextflow configuration <config-page>`.

(config-unscoped)=

## Unscoped options

`bucketDir`
: The remote work directory used by hybrid workflows. Equivalent to the `-bucket-dir` option of the `run` command.

`cleanup`
: If `true`, on a successful completion of a run all files in *work* directory are automatically deleted.

  :::{warning}
  The use of the `cleanup` option will prevent the use of the *resume* feature on subsequent executions of that pipeline run.
  :::
  :::{warning}
  The `cleanup` option is not supported for remote work directories, such as Amazon S3, Google Cloud Storage, and Azure Blob Storage.
  :::

`dumpHashes`
: If `true`, dump task hash keys in the log file, for debugging purposes. Equivalent to the `-dump-hashes` option of the `run` command.

`outputDir`
: :::{versionadded} 24.10.0
  :::
: Defines the pipeline output directory. Equivalent to the `-output-dir` option of the `run` command.

`resume`
: If `true`, enable the use of previously cached task executions. Equivalent to the `-resume` option of the `run` command.

`workDir`
: Defines the pipeline work directory. Equivalent to the `-work-dir` option of the `run` command.

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

Read the {ref}`container-apptainer` page to learn more about how to use Apptainer containers with Nextflow.

(config-aws)=

## `aws`

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
: :::{versionadded} 23.12.0-edge
  :::
: The AWS Batch Execution Role ARN that needs to be used to execute the Batch Job. This is mandatory when using AWS Fargate platform type. See [AWS documentation](https://docs.aws.amazon.com/batch/latest/userguide/execution-IAM-role.html) for more details.

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
: :::{versionadded} 23.12.0-edge
  :::
: Allow specifying the compute platform type used by AWS Batch, that can be either `ec2` or `fargate`. See AWS documentation to learn more about [AWS Fargate platform type](https://docs.aws.amazon.com/batch/latest/userguide/fargate.html) for AWS Batch.

`aws.batch.retryMode`
: The retry mode configuration setting, to accommodate rate-limiting on [AWS services](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html) (default: `standard`, other options: `legacy`, `adaptive`); this handling is delegated to AWS. To have Nextflow handle retries instead, use `built-in`.

`aws.batch.schedulingPriority`
: :::{versionadded} 23.01.0-edge
  :::
: The scheduling priority for all tasks when using [fair-share scheduling for AWS Batch](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/) (default: `0`)

`aws.batch.shareIdentifier`
: :::{versionadded} 22.09.0-edge
  :::
: The share identifier for all tasks when using [fair-share scheduling for AWS Batch](https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/)

`aws.batch.terminateUnschedulableJobs`
: :::{versionadded} 25.03.0-edge
  :::
: When `true`, jobs that cannot be scheduled for lack of resources or misconfiguration are terminated automatically (default: `false`). The pipeline may complete with an error status depending on the error strategy defined for the corresponding jobs.

`aws.batch.volumes`
: One or more container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. `/host/path:/mount/path[:ro|rw]`. Multiple mounts can be specified separating them with a comma or using a list object.

`aws.client.anonymous`
: Allow the access of public S3 buckets without the need to provide AWS credentials (default: `false`). Any service that does not accept unsigned requests will return a service access error.

`aws.client.s3Acl`
: Allow the setting of predefined bucket permissions, also known as *canned ACL*. Permitted values are `Private`, `PublicRead`, `PublicReadWrite`, `AuthenticatedRead`, `LogDeliveryWrite`, `BucketOwnerRead`, `BucketOwnerFullControl`, and `AwsExecRead` (default: none). See [Amazon docs](https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl) for details.

`aws.client.connectionTimeout`
: The amount of time to wait (in milliseconds) when initially establishing a connection before timing out (default: `10000`).

`aws.client.endpoint`
: The AWS S3 API entry point e.g. `https://s3-us-west-1.amazonaws.com`. The endpoint must include the protocol prefix e.g. `https://`.

`aws.client.glacierAutoRetrieval`
: :::{deprecated} 24.02.0-edge
  Glacier auto-retrieval is no longer supported. Instead, consider using the AWS CLI to restore any Glacier objects before or at the beginning of your pipeline (i.e. in a Nextflow process).
  :::
: Enable auto retrieval of S3 objects with a Glacier storage class (default: `false`).

`aws.client.glacierExpirationDays`
: :::{deprecated} 24.02.0-edge
  :::
: The time, in days, between when an object is restored to the bucket and when it expires (default: `7`).

`aws.client.glacierRetrievalTier`
: :::{deprecated} 24.02.0-edge
  :::
: The retrieval tier to use when restoring objects from Glacier, one of [`Expedited`, `Standard`, `Bulk`].

`aws.client.maxConcurrency`
: :::{versionadded} 25.06.0-edge
  :::
: The maximum number of concurrent S3 transfers used by the S3 transfer manager. By default, this setting is determined by `aws.client.targetThroughputInGbps`. Modifying this value can affect the amount of memory used for S3 transfers.

`aws.client.maxConnections`
: The maximum number of open HTTP connections used by the S3 transfer manager (default: `50`).

`aws.client.maxErrorRetry`
: The maximum number of retry attempts for failed retryable requests (default: `-1`).

`aws.client.maxNativeMemory`
: :::{versionadded} 25.06.0-edge
  :::
: The maximum native memory used by the S3 transfer manager. By default, this setting is determined by `aws.client.targetThroughputInGbps`.

`aws.client.minimumPartSize`
: :::{versionadded} 25.06.0-edge
  :::
: The minimum part size used by the S3 transfer manager for multi-part uploads (default: `8 MB`).

`aws.client.multipartThreshold`
: :::{versionadded} 25.06.0-edge
  :::
: The object size threshold used by the S3 transfer manager for performing multi-part uploads (default: same as `aws.cllient.minimumPartSize`).

`aws.client.protocol`
: :::{deprecated} 25.06.0-edge
  This option is no longer supported.
  :::
: The protocol to use when connecting to AWS. Can be `http` or `https` (default: `'https'`).

`aws.client.proxyHost`
: The proxy host to connect through.

`aws.client.proxyPort`
: The port to use when connecting through a proxy.

`aws.client.proxyScheme`
: :::{versionadded} 25.06.0-edge
  :::
: The protocol scheme to use when connecting through a proxy. Can be `http` or `https` (default: `'http'`).

`aws.client.proxyUsername`
: The user name to use when connecting through a proxy.

`aws.client.proxyPassword`
: The password to use when connecting through a proxy.

`aws.client.requesterPays`
: :::{versionadded} 24.05.0-edge
  :::
: Use [Requester Pays](https://docs.aws.amazon.com/AmazonS3/latest/userguide/RequesterPaysBuckets.html) for S3 buckets (default: `false`).

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

`aws.client.storageEncryption`
: The S3 server side encryption to be used when saving objects on S3. Can be `AES256` or `aws:kms` (default: none).

`aws.client.storageKmsKeyId`
: :::{versionadded} 22.05.0-edge
  :::
: The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket.

`aws.client.targetThroughputInGbps`
: :::{versionadded} 25.06.0-edge
  :::
: The target network throughput (in Gbps) used by the S3 transfer manager (default: `10`). This setting is not used when `aws.client.maxConcurrency` and `aws.client.maxNativeMemory` are specified.

`aws.client.transferManagerThreads`
: :::{versionadded} 25.06.0-edge
  :::
: Number of threads used by the S3 transfer manager (default `10`).

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

Read the {ref}`azure-page` page for more information.

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

`azure.batch.accountName`
: The batch service account name. Defaults to environment variable `AZURE_BATCH_ACCOUNT_NAME`.

`azure.batch.accountKey`
: The batch service account key. Defaults to environment variable `AZURE_BATCH_ACCOUNT_KEY`.

`azure.batch.allowPoolCreation`
: Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: `false`).

`azure.batch.autoPoolMode`
: Enable the automatic creation of batch pools depending on the pipeline resources demand (default: `true`).

`azure.batch.copyToolInstallMode`
: Specify where the `azcopy` tool used by Nextflow. When `node` is specified it's copied once during the pool creation. When `task` is provider, it's installed for each task execution. Finally when `off` is specified, the `azcopy` tool is not installed (default: `node`).

`azure.batch.deleteJobsOnCompletion`
: Delete all jobs when the workflow completes (default: `false`).
: :::{versionchanged} 23.08.0-edge
  Default value was changed from `true` to `false`.
  :::

`azure.batch.deletePoolsOnCompletion`
: Delete all compute node pools when the workflow completes (default: `false`).

`azure.batch.deleteTasksOnCompletion`
: :::{versionadded} 23.08.0-edge
  :::
: Delete each task when it completes (default: `true`).
: Although this setting is enabled by default, failed tasks will not be deleted unless it is explicitly enabled. This way, the default behavior is that successful tasks are deleted while failed tasks are preserved for debugging purposes.

`azure.batch.endpoint`
: The batch service endpoint e.g. `https://nfbatch1.westeurope.batch.azure.com`.

`azure.batch.location`
: The name of the batch service region, e.g. `westeurope` or `eastus2`. This is not needed when the endpoint is specified.

`azure.batch.jobMaxWallClockTime`
: :::{versionadded} 25.04.0
  :::
: The maximum elapsed time that jobs may run, measured from the time they are created. If jobs do not complete within this time limit, the Batch service terminates them and any tasks still running (default: `30d`).

`azure.batch.terminateJobsOnCompletion`
: :::{versionadded} 23.05.0-edge
  :::
: When the workflow completes, set all jobs to terminate on task completion. (default: `true`).

`azure.batch.pools.<name>.autoScale`
: Enable autoscaling feature for the pool identified with `<name>`.

`azure.batch.pools.<name>.fileShareRootPath`
: If mounting File Shares, this is the internal root mounting point. Must be `/mnt/resource/batch/tasks/fsmounts` for CentOS nodes or `/mnt/batch/tasks/fsmounts` for Ubuntu nodes (default is for CentOS).

`azure.batch.pools.<name>.lowPriority`
: Enable the use of low-priority VMs (default: `false`).

`azure.batch.pools.<name>.maxVmCount`
: Specify the max of virtual machine when using auto scale option.

`azure.batch.pools.<name>.mountOptions`
: Specify the mount options for mounting the file shares (default: `-o vers=3.0,dir_mode=0777,file_mode=0777,sec=ntlmssp`).

`azure.batch.pools.<name>.offer`
: Specify the offer type of the virtual machine type used by the pool identified with `<name>` (default: `centos-container`).

`azure.batch.pools.<name>.privileged`
: Enable the task to run with elevated access. Ignored if `runAs` is set (default: `false`).

`azure.batch.pools.<name>.publisher`
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
: Specify the ID of the Compute Node agent SKU which the pool identified with `<name>` supports (default: `batch.node.centos 8`).

`azure.batch.pools.<name>.startTask.script`
: :::{versionadded} 24.03.0-edge
  :::
: Specify the `startTask` that is executed as the node joins the Azure Batch node pool.

`azure.batch.pools.<name>.startTask.privileged`
: :::{versionadded} 24.03.0-edge
  :::
: Enable the `startTask` to run with elevated access (default: `false`).

`azure.batch.pools.<name>.virtualNetwork`
: :::{versionadded} 23.03.0-edge
  :::
: Specify the subnet ID of a virtual network in which to create the pool.

`azure.batch.pools.<name>.vmCount`
: Specify the number of virtual machines provisioned by the pool identified with `<name>`.

`azure.batch.pools.<name>.vmType`
: Specify the virtual machine type used by the pool identified with `<name>`.

`azure.batch.poolIdentityClientId`
: :::{versionadded} 25.05.0-edge
  :::
: Specify the client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) that is available on all Azure Batch node pools. This identity will be used by Fusion to authenticate to Azure storage. If set to `'auto'`, Fusion will use the first available managed identity.

`azure.managedIdentity.clientId`
: Specify the client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview). See {ref}`azure-managed-identities` for more details. Defaults to environment variable `AZURE_MANAGED_IDENTITY_USER`.

`azure.managedIdentity.system`
: When `true`, uses the system-assigned [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) to authenticate Azure resources. See {ref}`azure-managed-identities` for more details. Defaults to environment variable `AZURE_MANAGED_IDENTITY_SYSTEM`.

`azure.registry.server`
: Specify the container registry from which to pull the Docker images (default: `docker.io`).

`azure.registry.userName`
: Specify the username to connect to a private container registry.

`azure.registry.password`
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
: The blob storage account name. Defaults to environment variable `AZURE_STORAGE_ACCOUNT_NAME`.

`azure.storage.accountKey`
: The blob storage account key. Defaults to environment variable `AZURE_STORAGE_ACCOUNT_KEY`.

`azure.storage.sasToken`
: The blob storage shared access signature token, which can be provided as an alternative to `accountKey`. Defaults to environment variable `AZURE_STORAGE_SAS_TOKEN`.

`azure.storage.tokenDuration`
: The duration of the shared access signature token created by Nextflow when the `sasToken` option is *not* specified (default: `48h`).

(config-charliecloud)=

## `charliecloud`

The `charliecloud` scope controls how [Charliecloud](https://hpc.github.io/charliecloud/) containers are executed by Nextflow.

The following settings are available:

`charliecloud.enabled`
: Execute tasks with Charliecloud containers (default: `false`).

`charliecloud.writeFake`
: Enable `writeFake` with charliecloud (default: `true`) This allows to run containers from storage in writeable mode, using overlayfs. `writeFake` requires unprivileged `overlayfs` (Linux kernel >= 5.11). For full support, tempfs with xattrs in the user namespace (Linux kernel >= 6.6) is required, see [charliecloud documentation](https://hpc.github.io/charliecloud/ch-run.html#ch-run-overlay) for details.

`charliecloud.cacheDir`
: The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

`charliecloud.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`charliecloud.pullTimeout`
: The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: `20 min`).

`charliecloud.runOptions`
: Specify extra command line options supported by the `ch-run` command.

`charliecloud.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

`charliecloud.registry`
: The registry from where images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

Read the {ref}`container-charliecloud` page to learn more about how to use Charliecloud containers with Nextflow.

(config-conda)=

## `conda`

The `conda` scope controls the creation of Conda environments by the Conda package manager.

The following settings are available:

`conda.enabled`
: Enables the use of Conda environments (default: `false`).

`conda.cacheDir`
: Defines the path where Conda environments are stored. Ensure the path is accessible from all compute nodes when using a shared file system.

`conda.channels`
: Defines the Conda channels that can be used to resolve Conda packages. Channels can be defined as a list (e.g., `['bioconda','conda-forge']`) or a comma separated list string (e.g., `'bioconda,conda-forge'`). Channel priority decreases from left to right.

`conda.createOptions`
: Defines extra command line options supported by the `conda create` command. See the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/commands/create.html) for more information.

`conda.createTimeout`
: Defines the amount of time the Conda environment creation can last (default: `20 min`). The creation process is terminated when the timeout is exceeded.

`conda.useMamba`
: Uses the `mamba` binary instead of `conda` to create the Conda environments (default: `false`). See the [Mamba documentation](https://github.com/mamba-org/mamba) for more information about Mamba.

`conda.useMicromamba`
: :::{versionadded} 22.05.0-edge
  :::
: Uses the `micromamba` binary instead of `conda` to create Conda environments (default: `false`). See the [Micromamba documentation](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for more information about Micromamba.

See {ref}`conda-page` for more information about using Conda environments with Nextflow.

(config-dag)=

## `dag`

The `dag` scope controls the workflow diagram generated by Nextflow.

The following settings are available:

`dag.enabled`
: When `true` enables the generation of the DAG file (default: `false`).

`dag.depth`
: :::{versionadded} 23.10.0
  :::
: *Only supported by the HTML and Mermaid renderers.*
: Controls the maximum depth at which to render sub-workflows (default: no limit).

`dag.direction`
: :::{versionadded} 23.10.0
  :::
: *Supported by Graphviz, DOT, HTML and Mermaid renderers.*
: Controls the direction of the DAG, can be `'LR'` (left-to-right) or `'TB'` (top-to-bottom) (default: `'TB'`).

`dag.file`
: Graph file name (default: `dag-<timestamp>.html`).

`dag.overwrite`
: When `true` overwrites any existing DAG file with the same name (default: `false`).

`dag.verbose`
: :::{versionadded} 23.10.0
  :::
: *Only supported by the HTML and Mermaid renderers.*
: When `false`, channel names are omitted, operators are collapsed, and empty workflow inputs are removed (default: `false`).

Read the {ref}`workflow-diagram` page to learn more about the workflow graph that can be generated by Nextflow.

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
: Fix ownership of files created by the docker container (default: `false`).

`docker.legacy`
: Use command line options removed since Docker 1.10.0 (default: `false`).

`docker.mountFlags`
: Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`.

`docker.registry`
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`docker.registryOverride`
: :::{versionadded} 25.06.0-edge
  :::
: When `true`, forces the override of the registry name in fully qualified container image names with the registry specified by `docker.registry` (default: `false`). This setting allows you to redirect container image pulls from their original registry to a different registry, such as a private mirror or proxy.

`docker.remove`
: Clean-up the container after the execution (default: `true`). See the [Docker documentation](https://docs.docker.com/engine/reference/run/#clean-up---rm) for details.

`docker.runOptions`
: Specify extra command line options supported by the `docker run` command. See the [Docker documentation](https://docs.docker.com/engine/reference/run/) for details.

`docker.sudo`
: Executes Docker run command as `sudo` (default: `false`).

`docker.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

`docker.tty`
: Allocates a pseudo-tty (default: `false`).

Read the {ref}`container-docker` page to learn more about how to use Docker containers with Nextflow.

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
: The maximum number of CPUs made available by the underlying system. Used only by the `local` executor.

`executor.dumpInterval`
: Determines how often to log the executor status (default: `5min`).

`executor.exitReadTimeout`
: *Used only by grid executors.*
: Determines how long to wait before returning an error status when a process is terminated but the `.exitcode` file does not exist or is empty (default: `270 sec`).

`executor.jobName`
: Determines the name of jobs submitted to the underlying cluster executor e.g. `executor.jobName = { "$task.name - $task.hash" }`. Make sure the resulting job name matches the validation constraints of the underlying batch scheduler.
: This setting is supported by the following executors: Bridge, Condor, Flux, HyperQueue, Lsf, Moab, Nqsii, Oar, PBS, PBS Pro, SGE, SLURM and Google Batch.

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
: *Used only by the {ref}`lsf-executor` executor.*
: Specifies Platform LSF *per-job* memory limit mode (default: `false`).

`executor.perTaskReserve`
: *Used only by the {ref}`lsf-executor` executor.*
: Specifies Platform LSF *per-task* memory reserve mode (default: `false`).

`executor.pollInterval`
: Defines the polling frequency for process termination detection. Default varies for each executor (see below).

`executor.queueGlobalStatus`
: :::{versionadded} 23.01.0-edge
  :::
: Determines how job status is retrieved. When `false` only the queue associated with the job execution is queried. When `true` the job status is queried globally i.e. irrespective of the submission queue (default: `false`).

`executor.queueSize`
: The number of tasks the executor will handle in a parallel manner. A queue size of zero corresponds to no limit. Default varies for each executor (see below).

`executor.queueStatInterval`
: *Used only by grid executors.*
: Determines how often to fetch the queue status from the scheduler (default: `1min`).

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

`executor.retry.maxAttempt`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Max attempts when retrying failed job submissions (default: `3`).

`executor.retry.maxDelay`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Max delay when retrying failed job submissions (default: `30s`).

`executor.submit.retry.reason`
: :::{versionadded} 22.03.0-edge
  :::
: *Used only by grid executors.*
: Regex pattern that when verified cause a failed submit operation to be re-tried (default: `Socket timed out`).

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

`fusion.enabled`
: Enable Fusion file system (default: `false`).

`fusion.cacheSize`
: :::{versionadded} 23.11.0-edge
  :::
: Fusion client local cache size limit.

`fusion.containerConfigUrl`
: URL for downloading the container layer provisioning the Fusion client. 

`fusion.exportStorageCredentials`
: :::{versionadded} 23.05.0-edge
  Previously named `fusion.exportAwsAccessKeys`.
  :::
: Enable access to credentials for the underlying object storage are exported to the task environment (default: `false`).

`fusion.logLevel`
: Fusion client log level.

`fusion.logOutput`
: Log output location.

`fusion.privileged`
: :::{versionadded} 23.10.0
  :::
: Enable privileged containers for Fusion (default: `true`)
: Non-privileged use is supported only on Kubernetes with the [k8s-fuse-plugin](https://github.com/nextflow-io/k8s-fuse-plugin) or a similar FUSE device plugin.

`fusion.snapshots`
: :::{versionadded} 25.03.0-edge
  :::
: *Currently only supported for AWS Batch*
: Enable Fusion snapshotting (preview, default: `false`). This feature allows Fusion to automatically restore a job when it is interrupted by a spot reclamation.

`fusion.tags`
: Pattern for applying tags to files created via the Fusion client (default: `[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)`).
: Set to `false` to disable.

(config-google)=

## `google`

The `google` scope allows you to configure the interactions with Google Cloud, including Google Cloud Batch and Google Cloud Storage.

Read the {ref}`google-page` page for more information.

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
: The Google Cloud project ID to use for pipeline execution.

`google.batch.allowedLocations`
: :::{versionadded} 22.12.0-edge
  :::
: Define the set of allowed locations for VMs to be provisioned. See [Google documentation](https://cloud.google.com/batch/docs/reference/rest/v1/projects.locations.jobs#locationpolicy) for details (default: no restriction).

`google.batch.autoRetryExitCodes`
: :::{versionadded} 24.07.0-edge
  :::
: Defines the list of exit codes that will trigger Google Batch to automatically retry the job (default: `[50001]`). For this setting to take effect, `google.batch.maxSpotAttempts` must be greater than 0. See [Google Batch documentation](https://cloud.google.com/batch/docs/troubleshooting#reserved-exit-codes) for the complete list of retryable exit codes.

`google.batch.bootDiskImage`
: :::{versionadded} 24.08.0-edge
  :::
: Set the image URI of the virtual machine boot disk, e.g `batch-debian`. See [Google documentation](https://cloud.google.com/batch/docs/vm-os-environment-overview#vm-os-image-options) for details (default: none).

`google.batch.bootDiskSize`
: Set the size of the virtual machine boot disk, e.g `50.GB` (default: none).

`google.batch.cpuPlatform`
: Set the minimum CPU Platform, e.g. `'Intel Skylake'`. See [Specifying a minimum CPU Platform for VM instances](https://cloud.google.com/compute/docs/instances/specify-min-cpu-platform#specifications) (default: none).

`google.batch.gcsfuseOptions`
: :::{versionadded} 25.03.0-edge
  :::
: Defines a list of custom mount options for `gcsfuse` (default: `['-o rw', '-implicit-dirs']`).

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

  You can specify the network as a full or partial URL. For example, the following are all valid URLs:

  - https://www.googleapis.com/compute/v1/projects/{project}/global/networks/{network}
  - projects/{project}/global/networks/{network}
  - global/networks/{network}

`google.batch.networkTags`
: The network tags to be applied to the instances created by Google Batch jobs. Network tags are used to apply firewall rules and control network access (e.g., `['allow-ssh', 'allow-http']`). 

: Network tags are ignored when using instance templates. See [Add network tags](https://cloud.google.com/vpc/docs/add-remove-network-tags) for more information.

`google.batch.serviceAccountEmail`
: Define the Google service account email to use for the pipeline execution. If not specified, the default Compute Engine service account for the project will be used.

: Note that the `google.batch.serviceAccountEmail` service account will only be used for spawned jobs, not for the Nextflow process itself.  See the [Google Cloud](https://www.nextflow.io/docs/latest/google.html#credentials) documentation for more information on credentials.

`google.batch.spot`
: When `true` enables the usage of *spot* virtual machines or `false` otherwise (default: `false`).

`google.batch.subnetwork`
: The URL of an existing subnetwork resource in the network to which the VM will be attached.

  You can specify the subnetwork as a full or partial URL. For example, the following are all valid URLs:

  - https://www.googleapis.com/compute/v1/projects/{project}/regions/{region}/subnetworks/{subnetwork}
  - projects/{project}/regions/{region}/subnetworks/{subnetwork}
  - regions/{region}/subnetworks/{subnetwork}

`google.batch.usePrivateAddress`
: When `true` the VM will NOT be provided with a public IP address, and only contain an internal IP. If this option is enabled, the associated job can only load docker images from Google Container Registry, and the job executable cannot use external services other than Google APIs (default: `false`).

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
: Automatically mounts host paths into the task pods (default: `false`). Only intended for development purposes when using a single node.

`k8s.computeResourceType`
: :::{versionadded} 22.05.0-edge
  :::
: Define whether use Kubernetes `Pod` or `Job` resource type to carry out Nextflow tasks (default: `Pod`).

`k8s.context`
: Defines the Kubernetes [configuration context name](https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/) to use.

`k8s.cpuLimits`
: :::{versionadded} 24.04.0
  :::
: When `true`, set both the pod CPUs `request` and `limit` to the value specified by the `cpus` directive, otherwise set only the `request` (default: `false`).
: This setting is useful when a K8s cluster requires a CPU limit to be defined through a [LimitRange](https://kubernetes.io/docs/concepts/policy/limit-range/).

`k8s.debug.yaml`
: When `true`, saves the pod spec for each task to `.command.yaml` in the task directory (default: `false`).

`k8s.fetchNodeName`
: :::{versionadded} 22.05.0-edge
  :::
: If you trace the hostname, activate this option (default: `false`).

`k8s.fuseDevicePlugin`
: :::{versionadded} 24.01.0-edge
  :::
: The FUSE device plugin to be used when enabling Fusion in unprivileged mode (default: `['nextflow.io/fuse': 1]`).

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
: Defines the Kubernetes API max request retries (default: `4`).

`k8s.namespace`
: Defines the Kubernetes namespace to use (default: `default`).

`k8s.pod`
: Allows the definition of one or more pod configuration options such as environment variables, config maps, secrets, etc. It allows the same settings as the {ref}`process-pod` process directive.
: When using the `kuberun` command, this setting also applies to the submitter pod.

`k8s.projectDir`
: Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: `<volume-claim-mount-path>/projects`).

`k8s.pullPolicy`
: Defines the strategy to be used to pull the container image e.g. `'Always'`.

`k8s.retryPolicy.delay`
: Delay when retrying failed API requests (default: `500ms`).

`k8s.retryPolicy.jitter`
: Jitter value when retrying failed API requests (default: `0.25`).

`k8s.retryPolicy.maxAttempts`
: Max attempts when retrying failed API requests (default: `4`).

`k8s.retryPolicy.maxDelay`
: Max delay when retrying failed API requests (default: `90s`).

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

(config-lineage)=

## `lineage`

The `lineage` scope controls the generation of {ref}`cli-lineage` metadata.

The following settings are available:

`lineage.enabled`
: Enable generation of lineage metadata (default: `false`).

`lineage.store.location`
: Defines the location of the lineage metadata store (default: `./.lineage`).

(config-mail)=

## `mail`

The `mail` scope controls the mail server used to send email notifications.

The following settings are available:

`mail.debug`
: Enables Java Mail logging for debugging purposes (default: `false`).

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

## `manifest`

The `manifest` scope allows you to define some meta-data information needed when publishing or running your pipeline.

The following settings are available:

`manifest.author`
: :::{deprecated} 24.09.0-edge
  Use `manifest.contributors` instead.
  :::
: Project author name (use a comma to separate multiple names).

`manifest.contributors`
: :::{versionadded} 24.09.0-edge
  :::
: List of project contributors. Should be a list of maps. The following fields are supported in the contributor map:
  - `name`: the contributor's name 
  - `affiliation`: the contributor's affiliated organization
  - `email`: the contributor's email address
  - `github`: the contributor's GitHub URL
  - `contribution`: list of contribution types, each element can be one of `'author'`, `'maintainer'`, or `'contributor'`
  - `orcid`: the contributor's [ORCID](https://orcid.org/) URL

`manifest.defaultBranch`
: Git repository default branch (default: `master`).

`manifest.description`
: Free text describing the workflow project.

`manifest.docsUrl`
: Project documentation URL.

`manifest.doi`
: Project related publication DOI identifier.

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

  This setting may be useful to ensure that a specific version is used:

  ```groovy
  manifest.nextflowVersion = '1.2.3'        // exact match
  manifest.nextflowVersion = '1.2+'         // 1.2 or later (excluding 2 and later)
  manifest.nextflowVersion = '>=1.2'        // 1.2 or later
  manifest.nextflowVersion = '>=1.2, <=1.5' // any version in the 1.2 .. 1.5 range
  manifest.nextflowVersion = '!>=1.2'       // with ! prefix, stop execution if current version does not match required version.
  ```

`manifest.organization`
: Project organization

`manifest.recurseSubmodules`
: Pull submodules recursively from the Git repository.

`manifest.version`
: Project version number.

Read the {ref}`sharing-page` page to learn how to publish your pipeline to GitHub, BitBucket or GitLab.

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

The `notification` scope allows you to define the automatic sending of a notification email message when the workflow execution terminates.

`notification.attributes`
: A map object modelling the variables that can be used in the template file.

`notification.enabled`
: Send a notification message when the workflow execution completes (default: `false`).

`notification.from`
: Sender address for the notification email message.

`notification.template`
: Path of a template file which provides the content of the notification message.

`notification.to`
: Recipient address for the notification email. Multiple addresses can be specified separating them with a comma.

The notification message is sent my using the STMP server defined in the configuration {ref}`mail scope<config-mail>`.

If no mail configuration is provided, it tries to send the notification message by using the external mail command eventually provided by the underlying system (e.g. `sendmail` or `mail`).

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
: Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`.

`podman.registry`
: The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`podman.remove`
: Clean-up the container after the execution (default: `true`).

`podman.runOptions`
: Specify extra command line options supported by the `podman run` command.

`podman.temp`
: Mounts a path of your choice as the `/tmp` directory in the container. Use the special value `auto` to create a temporary directory each time a container is created.

Read the {ref}`container-podman` page to learn more about how to use Podman containers with Nextflow.

(config-report)=

## `report`

The `report` scope allows you to configure the workflow {ref}`execution-report`.

The following settings are available:

`report.enabled`
: Create the execution report on workflow completion (default: `false`).

`report.file`
: The path of the created execution report file (default: `report-<timestamp>.html`).

`report.overwrite`
: When `true` overwrites any existing report file with the same name (default: `false`).

(config-sarus)=

## `sarus`

The ``sarus`` scope controls how [Sarus](https://sarus.readthedocs.io) containers are executed by Nextflow.

The following settings are available:

`sarus.enabled`
: Execute tasks with Sarus containers (default: `false`).

`sarus.envWhitelist`
: Comma separated list of environment variable names to be included in the container environment.

`sarus.runOptions`
: Specify extra command line options supported by the `sarus run` command. For details see the [Sarus user guide](https://sarus.readthedocs.io/en/stable/user/user_guide.html).

`sarus.tty`
: Allocates a pseudo-tty (default: `false`).

Read the {ref}`container-sarus` page to learn more about how to use Sarus containers with Nextflow.

(config-shifter)=

## `shifter`

The `shifter` scope controls how [Shifter](https://docs.nersc.gov/programming/shifter/overview/) containers are executed by Nextflow.

The following settings are available:

`shifter.enabled`
: Execute tasks with Shifter containers (default: `false`).

Read the {ref}`container-shifter` page to learn more about how to use Shifter containers with Nextflow.

(config-singularity)=

## `singularity`

The `singularity` scope controls how [Singularity](https://sylabs.io/singularity/) containers are executed by Nextflow.

The following settings are available:

`singularity.autoMounts`
: Automatically mounts host paths in the executed container (default: `true`). It requires the `user bind control` feature to be enabled in your Singularity installation.
: :::{versionchanged} 23.09.0-edge
  Default value was changed from `false` to `true`.
  :::

`singularity.cacheDir`
: The directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.

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
: When enabled, OCI (and Docker) container images are pull and converted to a SIF image file format implicitly by the Singularity run command, instead of Nextflow. Requires Singularity 3.11 or later (default: `false`).

  :::{note}
  Leave `ociAutoPull` disabled if willing to build a Singularity native image with Wave (see the {ref}`wave-singularity` section).
  :::

`singularity.ociMode`
: :::{versionadded} 23.12.0-edge
  :::
: Enable OCI-mode, that allows running native OCI compliant container image with Singularity using `crun` or `runc` as low-level runtime. Note: it requires Singularity 4 or later. See `--oci` flag in the [Singularity documentation](https://docs.sylabs.io/guides/4.0/user-guide/oci_runtime.html#oci-mode) for more details and requirements (default: `false`).

  :::{note}
  Leave `ociMode` disabled if you are willing to build a Singularity native image with Wave (see the {ref}`wave-singularity` section).
  :::

`singularity.pullTimeout`
: The amount of time the Singularity pull can last, exceeding which the process is terminated (default: `20 min`).

`singularity.registry`
: :::{versionadded} 22.12.0-edge
  :::
: The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. `http://`.

`singularity.runOptions`
: Specify extra command line options supported by `singularity exec`.

Read the {ref}`container-singularity` page to learn more about how to use Singularity containers with Nextflow.

(config-spack)=

## `spack`

The `spack` scope controls the creation of a Spack environment by the Spack package manager.

The following settings are available:

`spack.cacheDir`
: Defines the path where Spack environments are stored. When using a compute cluster make sure to provide a shared file system path accessible from all compute nodes.

`spack.checksum`
: Enables checksum verification for source tarballs (default: `true`). Only disable when requesting a package version not yet encoded in the corresponding Spack recipe.

`spack.createTimeout`
: Defines the amount of time the Spack environment creation can last (default: `60 min`). The creation process is terminated when the timeout is exceeded.

`spack.parallelBuilds`
: Sets number of parallel package builds (Spack default: coincides with number of available CPU cores).

Nextflow does not allow for fine-grained configuration of the Spack package manager. Instead, this has to be performed directly on the host Spack installation. For more information see the [Spack documentation](https://spack.readthedocs.io).

(config-timeline)=

## `timeline`

The `timeline` scope controls the execution timeline report generated by Nextflow.

The following settings are available:

`timeline.enabled`
: Create the timeline report on workflow completion file (default: `false`).

`timeline.file`
: Timeline file name (default: `timeline-<timestamp>.html`).

`timeline.overwrite`
: When `true` overwrites any existing timeline file with the same name (default: `false`).

(config-tower)=

## `tower`

The `tower` scope controls the settings for the [Seqera Platform](https://seqera.io) (formerly Tower Cloud).

The following settings are available:

`tower.accessToken`
: The unique access token specific to your account on an instance of Seqera Platform.
: Your `accessToken` can be obtained from your Seqera Platform instance in the [Tokens page](https://cloud.seqera.io/tokens).

`tower.enabled`
: Send workflow tracing and execution metrics to Seqera Platform (default: `false`).

`tower.endpoint`
: The endpoint of your Seqera Platform instance (default: `https://api.cloud.seqera.io`).

`tower.workspaceId`
: The ID of the Seqera Platform workspace where the run should be added (default: the launching user personal workspace).
: The workspace ID can also be specified using the environment variable `TOWER_WORKSPACE_ID` (config file has priority over the environment variable).

(config-trace)=

## `trace`

The `trace` scope controls the layout of the execution trace file generated by Nextflow.

The following settings are available:

`trace.enabled`
: Create the execution trace file on workflow completion (default: `false`).

`trace.fields`
: Comma separated list of fields to be included in the report. The available fields are listed at {ref}`this page <trace-fields>`.

`trace.file`
: Trace file name (default: `trace-<timestamp>.txt`).

`trace.overwrite`
: When `true` overwrites any existing trace file with the same name (default: `false`).

`trace.raw`
: When `true` turns on raw number report generation i.e. date and time are reported as milliseconds and memory as number of bytes (default: `false`).

`trace.sep`
: Character used to separate values in each row (default: `\t`).

Read the {ref}`trace-report` page to learn more about the execution report that can be generated by Nextflow.

(config-wave)=

## `wave`

The `wave` scope provides advanced configuration for the use of {ref}`Wave containers <wave-page>`.

The following settings are available:

`wave.enabled`
: Enable the use of Wave containers (default: `false`).

`wave.build.repository`
: The container repository where images built by Wave are uploaded (note: the corresponding credentials must be provided in your Seqera Platform account).

`wave.build.cacheRepository`
: The container repository used to cache image layers built by the Wave service (note: the corresponding credentials must be provided in your Seqera Platform account).

`wave.build.compression.mode`
: :::{versionadded} 25.05.0-edge
  :::
: Defines the compression algorithm that should be used when building the container. Allowed values are: `gzip`, `estargz` and `zstd` (default: `gzip`).

`wave.build.compression.level`
: :::{versionadded} 25.05.0-edge
  :::
: Level of compression used when building a container depending the chosen algorithm: gzip, estargz (0-9) and zstd (0-22).
   .
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

`wave.endpoint`
: The Wave service endpoint (default: `https://wave.seqera.io`).

`wave.freeze`
: :::{versionadded} 23.07.0-edge
  :::
: Enables Wave container freezing (default: `false`). Wave will provision a non-ephemeral container image that will be pushed to a container repository of your choice. It requires the use of the `wave.build.repository` setting.
: It is also recommended to specify a custom cache repository using `wave.build.cacheRepository`.
: :::{note}
  The container repository authentication must be managed by the underlying infrastructure.
  :::

`wave.httpClient.connectTime`
: :::{versionadded} 22.06.0-edge
  :::
: Sets the connection timeout duration for the HTTP client connecting to the Wave service (default: `30s`).

`wave.mirror`
: :::{versionadded} 24.09.1-edge
  :::
: Enables Wave container mirroring (default: `false`).
: This feature allow mirroring (i.e. copying) the containers defined in your pipeline
  configuration to a container registry of your choice, so that pipeline tasks will pull the copied containers from the
  target registry instead of the original one.
: The resulting copied containers will maintain the name, digest and metadata.
: The target registry is expected to be specified by using the `wave.build.repository` option.
: :::{note}
  * This feature is only compatible with `wave.strategy = 'container'` option.
  * This feature cannot be used with Wave *freeze* mode.
  * The authentication of the resulting container images must be managed by the underlying infrastructure.
  :::

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

`wave.scan.mode`
: :::{versionadded} 24.09.1-edge
  :::
: Determines the container security scanning execution modality.

: This feature allows scanning for security vulnerability the container used in your pipeline. The following options can be specified:

  * `none`: No security scan is performed on the containers used by your pipeline.
  * `async`: The containers used by your pipeline are scanned for security vulnerability. The task execution is carried out independently of the security scan result.
  * `required`: The containers used by your pipeline are scanned for security vulnerability. The task is only executed if the corresponding container is not affected by a security vulnerability.

`wave.scan.allowedLevels`
: :::{versionadded} 24.09.1-edge
  :::
: Determines the allowed security levels when scanning containers for security vulnerabilities.

: Allowed values are: `low`, `medium`, `high`, `critical`. For example: `wave.scan.allowedLevels = 'low,medium'`.

: This option requires the use of `wave.scan.mode = 'required'`.

`wave.strategy`
: The strategy to be used when resolving ambiguous Wave container requirements (default: `'container,dockerfile,conda'`).

(config-workflow)=

## `workflow`

:::{versionadded} 24.10.0
:::

The `workflow` scope provides workflow execution options.

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
: The file publishing method (default: `'symlink'`). The following options are available:

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
: When `true` any existing file in the specified folder will be overwritten (default: `'standard'`). The following options are available:

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
: Specify arbitrary tags for published files. For example:
  ```groovy
  tags FOO: 'hello', BAR: 'world'
  ```
