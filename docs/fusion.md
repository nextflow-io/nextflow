(fusion-page)=

# Fusion file system

:::{versionadded} 22.10.0
:::

:::{versionadded} 23.02.0-edge
Support for Google Cloud Storage.
:::

Cloud object stores such as AWS S3 are scalable and cost-effective, but they don't present a POSIX (Portable Operating System Interface). This means containerized applications must copy data to and from cloud storage for every task — a slow and inefficient process.

Fusion is a virtual, lightweight, distributed file system that bridges the gap between pipelines and cloud-native storage. Fusion enables seamless filesystem I/O to cloud object stores via a standard POSIX interface, resulting in simpler pipeline logic and faster, more efficient pipeline execution.

See [Fusion file system](https://docs.seqera.io/fusion) for more information about Fusion features.

:::{note}
Fusion requires a license for use in Seqera Platform compute environments or directly in Nextflow. Fusion can be trialed at no cost. [Contact Seqera](https://seqera.io/contact-us/) for more details.
:::

## Get started

Fusion integrates with Nextflow directly and does not require any installation or change in pipeline code. It only requires use of a container runtime or a container computing service, such as Kubernetes, AWS Batch, or Google Cloud Batch.

To enable Fusion in your Nextflow pipeline, add the following to your `nextflow.config` file:

```{code-block} nextflow
:class: copyable
fusion.enabled = true
wave.enabled = true
tower.accessToken = '<PLATFORM_ACCESS_TOKEN>'
```

Replace `<PLATFORM_ACCESS_TOKEN>` with your Platform access token.

See [Get started](https://docs.seqera.io/fusion/get-started) for more information and guides to get started with Fusion.
