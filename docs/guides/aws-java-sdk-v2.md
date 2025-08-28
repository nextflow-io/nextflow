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

## S3 concurrency

You can use the `aws.client.maxConnections` config option to control the maximum number of concurrent HTTP connections to S3.

You can also use the `aws.client.targetThroughputInGbps` option to control the concurrency of S3 uploads and downloads specifically, based on the available network bandwidth. This setting is `10` by default, which means that Nextflow performs S3 transfers concurrently up to 10 Gbps of network throughput, regardless of the connection limit. All other S3 API calls are controlled by the connection limit.

Use these settings with virtual threads to achieve optimal performance for your environment. Increasing these settings beyond their defaults may improve performance for large runs. You can enable virtual threads by setting the `NXF_ENABLE_VIRTUAL_THREADS` environment variable to `true`.

## Multi-part uploads

Nextflow uploads large files to S3 as multi-part uploads. You can use the `aws.client.minimumPartSize` and `aws.client.multipartThreshold` config options to control when and how multi-part uploads are performed.

The following multi-part upload config options are no longer supported:

- `aws.client.uploadChunkSize`
- `aws.client.uploadMaxAttempts`
- `aws.client.uploadMaxThreads`
- `aws.client.uploadRetrySleep`
