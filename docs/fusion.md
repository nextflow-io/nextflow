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

Use Fusion directly in Seqera Platform compute environments, or add Fusion to your Nextflow pipeline configuration. See [Get started](https://docs.seqera.io/fusion/get-started) for more information about configuring Fusion for your compute environment.

## Configuration options

Add Fusion configuration options to your `nextflow.config` file. See {ref}`config-fusion` for a full list of Fusion configuration options.
