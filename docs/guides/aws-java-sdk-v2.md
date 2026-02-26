(aws-java-sdk-v2-page)=

# AWS Java SDK v2

AWS Java SDK v1 is reaching end of life at the end of 2025. Starting in version `25.06.0-edge`, Nextflow uses AWS Java SDK v2 in the `nf-amazon` plugin.

This migration introduced several breaking changes to the `aws.client` config scope, including new options and removed options. This page describes these changes and how they affect your Nextflow configuraiton.

## New HTTP client

The HTTP client used by SDK v2 does not support overriding certain advanced HTTP options. As a result, the following config options are no longer supported:

- `aws.client.protocol`
- `aws.client.signerOverride`
- `aws.client.socketRecvBufferSizeHint`
- `aws.client.socketSendBufferSizeHint`
- `aws.client.userAgent`

## S3 transfer manager

The *S3 transfer manager* is a subsystem of SDK v2 which handles S3 transfers, including S3 uploads and downloads.

The concurrency and throughput of the S3 transfer manager can be configured manually using the `aws.client.maxConcurrency` and `aws.client.maxNativeMemory` config options. Alternatively, the `aws.client.targetThroughputInGbps` config option can be used to set the previous two options automatically based on a target throughput.

## Multi-part uplaods

Multi-part uploads are handled by the S3 transfer manager. The `aws.client.minimumPartSize` and `aws.client.multipartThreshold` config options can be used to control when and how multi-part uploads are performed.

The following multi-part upload options are no longer supported:

- `aws.client.uploadChunkSize`
- `aws.client.uploadMaxAttempts`
- `aws.client.uploadMaxThreads`
- `aws.client.uploadRetrySleep`
