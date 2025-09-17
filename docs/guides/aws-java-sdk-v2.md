(aws-java-sdk-v2-page)=

# AWS Java SDK v2

AWS Java SDK v1 will reach end of life at the end of 2025. Starting with version `25.06.0-edge`, Nextflow uses AWS Java SDK v2 in the `nf-amazon` plugin.

This migration introduces several breaking changes to the `aws.client` config scope, including new and removed options. This page describes these changes and how they affect your Nextflow configuration.

## New HTTP client

The HTTP client in SDK v2 does not support overriding certain advanced HTTP options. As a result, the following config options are no longer supported:

- `aws.client.protocol`
- `aws.client.signerOverride`
- `aws.client.socketRecvBufferSizeHint`
- `aws.client.socketSendBufferSizeHint`
- `aws.client.userAgent`

## S3 transfer manager

The *S3 transfer manager* is a subsystem of SDK v2 that handles S3 uploads and downloads.

You can configure the concurrency and throughput of the S3 transfer manager manually using the `aws.client.maxConcurrency` and `aws.client.maxNativeMemory` configuration options. Alternatively, you can use the `aws.client.targetThroughputInGbps` option to set both values automatically based on a target throughput.

## Multi-part transfers

Multi-part transfer are handled by the S3 transfer manager. You can use the `aws.client.minimumPartSize` and `aws.client.multipartThreshold` config options to control when and how multi-part transfers are performed. 
Concurrent multi-part downloads can consume large heap memory space due to the buffer size created per transfer. To avoid out of memory errors, the size consumed by these buffers is limited to 400 MB. You can use `aws.client.maxDownloadBuffer` to change this value.

The following multi-part upload config options are no longer supported:

- `aws.client.uploadChunkSize`
- `aws.client.uploadMaxAttempts`
- `aws.client.uploadMaxThreads`
- `aws.client.uploadRetrySleep`
