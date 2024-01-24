(fusion-page)=

# Fusion file system

:::{versionadded} 22.10.0
:::

:::{versionadded} 23.02.0-edge
Support for Google Cloud Storage.
:::

## Introduction

Fusion is a distributed virtual file system for cloud-native data pipeline and optimised for Nextflow workloads.

It bridges the gap between cloud-native storage and data analysis workflow by implementing a thin client that allows any existing application to access object storage using the standard POSIX interface, thus simplifying and speeding up most operations. Currently it supports AWS S3 and Google Cloud Storage.

## Getting started

### Requirements

Fusion file system is designed to work with containerised workloads, therefore it requires the use of a container engine such as Docker or a container native platform for the execution of your pipeline e.g. AWS Batch or Kubernetes. It also requires the use of {ref}`Wave containers<wave-page>`.

### AWS S3 configuration

The AWS S3 bucket should be configured with the following IAM permissions:

```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:ListBucket"
            ],
            "Resource": [
                "arn:aws:s3:::YOUR-BUCKET-NAME"
            ]
        },
        {
            "Action": [
                "s3:GetObject",
                "s3:PutObject",
                "s3:PutObjectTagging",
                "s3:DeleteObject"
            ],
            "Resource": [
                "arn:aws:s3:::YOUR-BUCKET-NAME/*"
            ],
            "Effect": "Allow"
        }
    ]
}
```

## Use cases

### Local execution with S3 bucket as work directory

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Nextflow local executor. This configuration requires the use of Docker (or similar container engine) for the execution of your pipeline tasks.

The AWS S3 bucket credentials should be made accessible via standard `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` environment variables.

The following configuration should be added in your Nextflow configuration file:

```groovy
docker {
    enabled = true
}

fusion {
    enabled = true
    exportStorageCredentials = true
}

wave {
    enabled = true
}
```

Then you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

Replace `<YOUR PIPELINE>` and `<YOUR BUCKET>` with a pipeline script and bucket or your choice, for example:

```bash
nextflow run https://github.com/nextflow-io/rnaseq-nf -work-dir s3://nextflow-ci/scratch
```

### AWS Batch execution with S3 bucket as work directory

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the AWS Batch executor. The use of Fusion makes obsolete the need to create and configure a custom AMI that includes the `aws` command line tool, when setting up the AWS Batch compute environment.

The configuration for this deployment scenario looks like the following:

```groovy
fusion {
    enabled = true
}

wave {
    enabled = true
}

process {
    executor = 'awsbatch'
    queue = '<YOUR BATCH QUEUE>'
}

aws {
    region = '<YOUR AWS REGION>'
}
```

Then you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

### Kubernetes execution with S3 bucket as work directory

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Kubernetes executor.

The use of Fusion makes obsolete the need to create and manage and separate persistent volume and shared file system in the Kubernetes cluster.

The configuration for this deployment scenario looks like the following:

```groovy
wave {
    enabled = true
}

fusion {
    enabled = true
}

process {
    executor = 'k8s'
}

k8s {
    context = '<YOUR K8S CONFIGURATION CONTEXT>'
    namespace = '<YOUR K8S NAMESPACE>'
    serviceAccount = '<YOUR K8S SERVICE ACCOUNT>'
}
```

The `k8s.context` represents the Kubernetes configuration context to be used for the pipeline execution. This setting can be omitted if Nextflow itself is run as a pod in the Kubernetes clusters.

The `k8s.namespace` represents the Kubernetes namespace where the jobs submitted by the pipeline execution should be executed.

The `k8s.serviceAccount` represents the Kubernetes service account that should be used to grant the execution permission to jobs launched by Nextflow. You can find more details how to configure it as the [following link](https://github.com/seqeralabs/wave-showcase/tree/master/example8).

Having the above configuration in place, you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

## NVMe storage

The Fusion file system implements a lazy download and upload algorithm that runs in the background to transfer files in parallel to and from object storage into a container-local temporary folder. This means that the performance of the temporary folder inside the container (`/tmp` in a default setup) is key to achieving maximum performance.

The temporary folder is used only as a temporary cache, so the size of the volume can be much lower than the actual needs of your pipeline processes. Fusion has a built-in garbage collector that constantly monitors remaining disk space on the temporary folder and immediately evicts old cached entries when necessary.

The recommended setup to get maximum performance is to mount a NVMe disk as the temporary folder and run the pipeline with the {ref}`scratch <process-scratch>` directive set to `false` to also avoid stage-out transfer time.

Example configuration for using AWS Batch with [NVMe disks](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ssd-instance-store.html) to maximize performance:

```groovy
aws.batch.volumes = '/path/to/ec2/nvme:/tmp'
process.scratch = false
```

## Advanced settings

The following configuration options are available:

`fusion.enabled`
: Enable/disable the use of Fusion file system.

`fusion.cacheSize`
: :::{versionadded} 23.11.0-edge
:::
: The maximum size of the local cache used by the Fusion client.

`fusion.containerConfigUrl`
: The URL from where the container layer provisioning the Fusion client is downloaded.

`fusion.exportStorageCredentials`
: :::{versionadded} 23.05.0-edge
  This option was previously named `fusion.exportAwsAccessKeys`.
  :::
: When `true` the access credentials required by the underlying object storage are exported the pipeline jobs execution environment.

`fusion.logLevel`
: The level of logging emitted by the Fusion client.

`fusion.logOutput`
: Where the logging output is written. 

`fusion.privileged`
: :::{versionadded} 23.10.0
  :::
: This allows disabling the privileged container execution when using the Fusion file system.
  The effective use of this setting depends on the target execution. Currently, it's only supported by the Kubernetes
  executor which requires the use the [k8s-fuse-plugin](https://github.com/nextflow-io/k8s-fuse-plugin) to be installed
  in the target cluster (default: `true`).

`fusion.tags`
: The pattern that determines how tags are applied to files created via the Fusion client. To disable tags
  set it to `false`. (default: `[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)`)
