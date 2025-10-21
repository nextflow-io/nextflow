(aws-java-sdk-v2-page)=

# AWS Java SDK v2

AWS Java SDK v1 will reach end of life at the end of 2025. Starting in Nextflow 25.10, Nextflow uses AWS Java SDK v2 in the `nf-amazon` plugin.

This migration introduces several breaking changes to the `aws.client` config scope, including new and removed options. This page describes these changes and how they affect your Nextflow configuration.

## New HTTP client

The HTTP client in SDK v2 does not support overriding certain advanced HTTP options. As a result, the following config options are no longer supported:

- `aws.client.protocol`
- `aws.client.signerOverride`
- `aws.client.socketRecvBufferSizeHint`
- `aws.client.socketSendBufferSizeHint`
- `aws.client.userAgent`

## Parallel S3 operations

Nextflow manages S3 transfers, including uploads, downloads, and S3-to-S3 copies, separately from other S3 API calls such as listing a directory or retrieving object metadata.

Use the `aws.client.targetThroughputInGbps` option to control the concurrency of S3 transfers based on the available network bandwidth. This setting defaults to `10`, which allows Nextflow to perform concurrent S3 transfers up to 10 Gbps of network throughput.

Use the `aws.client.maxConnections` config option to control the maximum number of concurrent HTTP connections for all other S3 API calls.

## Multi-part transfers

Nextflow transfers large files to and from S3 as multipart transfers. Use the `aws.client.minimumPartSize` and `aws.client.multipartThreshold` configuration options to control when and how multipart transfers are performed.

Concurrent multipart downloads can consume a large amount of heap memory due to the buffer allocated by each transfer. To avoid out-of-memory errors, the size consumed by these buffers is limited to 400 MB by default. Use the `aws.client.maxDownloadHeapMemory` option to control this value.

The following multipart upload config options are no longer supported:

- `aws.client.uploadChunkSize`
- `aws.client.uploadMaxAttempts`
- `aws.client.uploadMaxThreads`
- `aws.client.uploadRetrySleep`
