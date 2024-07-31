(fusion-page)=

# Fusion file system

:::{versionadded} 22.10.0
:::

:::{versionadded} 23.02.0-edge
Support for Google Cloud Storage.
:::

## Introduction

Fusion is a distributed virtual file system for cloud-native data pipeline and optimised for Nextflow workloads.

It bridges the gap between cloud-native storage and data analysis workflow by implementing a thin client that allows any existing application to access object storage using the standard POSIX interface, thus simplifying and speeding up most operations.
Currently it supports AWS S3, Google Cloud Storage and Azure Blob containers.

## Getting started

The Fusion file system implements a lazy download and upload algorithm that runs in the background to transfer files in
parallel to and from object storage into a container-local temporary folder. This means that the performance of the disk
volume used to carry out your computation is key to achieving maximum performance.

By default Fusion uses the container `/tmp` directory as a temporary cache, so the size of the volume can be much lower
than the actual needs of your pipeline processes. Fusion has a built-in garbage collector that constantly monitors remaining
disk space on the temporary folder and immediately evicts old cached entries when necessary.

### Requirements

Fusion file system is designed to work with containerised workloads, therefore it requires the use of a container engine
such as Docker or a container native platform for the execution of your pipeline e.g. AWS Batch or Kubernetes. It also requires
the use of {ref}`Wave containers<wave-page>`.

### Azure Cloud

Fusion provides built-in support for [Azure Blob Storage](https://azure.microsoft.com/en-us/products/storage/blobs/)
when running in Azure Cloud.

The support for Azure does not require any specific setting other then enabling Wave and Fusion in your Nextflow
configuration. For example:

```
fusion.enabled = true
wave.enabled = true
process.executor = 'azure-batch'
tower.accessToken = '<your platform access token>'
```

Then run your pipeline using the usual command:

```
nextflow run <your pipeline> -work-dir az://<your blob container>/scratch
```

Azure machines come with fast SSDs attached, therefore no additional storage configuration is required however it is
recommended to use the machine types with larger data disks attached, denoted by the suffix `d` after the core number
(e.g. `Standard_E32*d*_v5`). These will increase the throughput of Fusion and reduce the chance of overloading the machine.

### AWS Cloud

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the AWS Batch executor.
The use of Fusion makes obsolete the need to create and configure a custom AMI that includes the `aws` command
line tool, when setting up the AWS Batch compute environment.

The configuration for this deployment scenario looks like the following:

```groovy
fusion.enabled = true
wave.enabled = true
process.executor = 'awsbatch'
process.queue = '<YOUR BATCH QUEUE>'
aws.region = '<YOUR AWS REGION>'
tower.accessToken = '<your platform access token>'
```

Then you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

For best performance make sure to use instance types that provide a NVMe disk as [instance storage](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/InstanceStorage.html).
If you are creating the AWS Batch compute environment by yourselves, you will need to make sure the NVMe is properly formatted (see below).


#### NVMe storage

The recommended setup to get maximum performance is to mount a NVMe disk as the temporary folder and run the pipeline with
the {ref}`scratch <process-scratch>` directive set to `false` to also avoid stage-out transfer time.

Example configuration for using AWS Batch with [NVMe disks](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ssd-instance-store.html) to maximize performance:

```groovy
aws.batch.volumes = '/path/to/ec2/nvme:/tmp'
process.scratch = false
```

:::{tip}
Seqera Platform is able to automatically format and configure the NVMe instance storage by enabling
the option "Use Fast storage" when creating the Batch compute environment.
:::

#### AWS IAM permissions

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

### Google Cloud

Fusion provides built-in support for [Google Storage](https://cloud.google.com/storage?hl=en)
when running in Google Cloud.

The support for Google does not require any specific setting other then enabling Wave and Fusion in your Nextflow
configuration. For example:

```
fusion.enabled = true
wave.enabled = true
process.executor = 'google-batch'
tower.accessToken = '<your platform access token>'
```

Then run your pipeline using the usual command:

```
nextflow run <your pipeline> -work-dir gs://<your google bucket>/scratch
```

When using Fusion, if the `process.disk` is not set, Nextflow will attach a single local SSD disk to the machine. The size of this disk can be much lower than the actual needs of your pipeline processes because Fusion uses it only as a temporal cache. Fusion is also compatible with other types of `process.disk`, but better performance is achieved when using local SSD disks.

### Kubernetes

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Kubernetes executor.

The use of Fusion makes obsolete the need to create and manage and separate persistent volume and shared file system in the Kubernetes cluster.

The configuration for this deployment scenario looks like the following:

```groovy
fusion.enabled = true
wave.enabled = true
process.executor = 'k8s'
k8s.context = '<YOUR K8S CONFIGURATION CONTEXT>'
k8s.namespace = '<YOUR K8S NAMESPACE>'
k8s.serviceAccount = '<YOUR K8S SERVICE ACCOUNT>'
tower.accessToken = '<your platform access token>'
```

The `k8s.context` represents the Kubernetes configuration context to be used for the pipeline execution. This setting can be omitted if Nextflow itself is run as a pod in the Kubernetes clusters.

The `k8s.namespace` represents the Kubernetes namespace where the jobs submitted by the pipeline execution should be executed.

The `k8s.serviceAccount` represents the Kubernetes service account that should be used to grant the execution permission to jobs launched by Nextflow. You can find more details how to configure it as the [following link](https://github.com/seqeralabs/wave-showcase/tree/master/example8).

Having the above configuration in place, you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

:::{note}
You an also use Fusion and Kubernetes with Azure Blob Storage and Google Storage using the same deployment approach.
:::

### Local execution

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Nextflow local executor. This configuration requires the use of Docker (or similar container engine) for the execution of your pipeline tasks.

The AWS S3 bucket credentials should be made accessible via standard `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY` environment variables.

The following configuration should be added in your Nextflow configuration file:

```groovy
docker.enabled = true
fusion.enabled = true
fusion.exportStorageCredentials = true
wave.enabled = true
```

Then you can run your pipeline using the following command:

```bash
nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch
```

Replace `<YOUR PIPELINE>` and `<YOUR BUCKET>` with a pipeline script and bucket or your choice, for example:

```bash
nextflow run https://github.com/nextflow-io/rnaseq-nf -work-dir s3://nextflow-ci/scratch
```

:::{warning}
The option `fusion.exportStorageCredentials` leaks the AWS credentials on the task launcher script created by Nextflow.
This option should only be used for development purposes.
:::

## Advanced settings

Fusion advanced configuration settings are described in the {ref}`Fusion <config-fusion>` section on the Nextflow configuration page.
