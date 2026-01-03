/*
 * Copyright 2020-2022, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

package nextflow.cloud.aws.config

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.spec.ConfigOption
import nextflow.script.dsl.Description
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Model AWS S3 client config settings. It is applied to all buckets when there is no a specific bucket configuration.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsS3ClientConfig extends AwsS3CommonConfig {

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) when initially establishing a connection before timing out (default: `10000`).
    """)
    final Integer connectionTimeout

    final Boolean debug

    /**
     * Maximum number of concurrent transfers used by S3 transfer manager. By default,
     * it is determined automatically by `targetThroughputInGbps`.
     */
    @ConfigOption
    final Integer maxConcurrency

    @ConfigOption
    @Description("""
        The maximum number of open HTTP connections used by the S3 client (default: `50`).
    """)
    final Integer maxConnections

    @ConfigOption
    @Description("""
        The maximum size for the heap memory buffer used by concurrent downloads. It must be at least 10 times the `minimumPartSize` (default:`400 MB`).
    """)
    final MemoryUnit maxDownloadHeapMemory

    @ConfigOption
    @Description("""
        The maximum number of retry attempts for failed retryable requests (default: `-1`).
    """)
    final Integer maxErrorRetry

    /**
     * Maximum native memory used by S3 transfer manager. By default, it is
     * determined automatically by `targetThroughputInGbps`.
     */
    @ConfigOption
    final MemoryUnit maxNativeMemory

    @ConfigOption
    @Description("""
        The minimum part size used for multipart uploads to S3 (default: `8 MB`).
    """)
    final MemoryUnit minimumPartSize

    @ConfigOption
    @Description("""
        The object size threshold used for multipart uploads to S3 (default: same as `aws.cllient.minimumPartSize`).
    """)
    final MemoryUnit multipartThreshold

    @ConfigOption
    @Description("""
        The proxy host to connect through.
    """)
    final String proxyHost

    @ConfigOption
    @Description("""
        The port to use when connecting through a proxy.
    """)
    final Integer proxyPort

    @ConfigOption
    @Description("""
        The protocol scheme to use when connecting through a proxy. Can be `http` or `https` (default: `'http'`).
    """)
    final String proxyScheme

    @ConfigOption
    @Description("""
        The user name to use when connecting through a proxy.
    """)
    final String proxyUsername

    @ConfigOption
    @Description("""
        The password to use when connecting through a proxy.
    """)
    final String proxyPassword

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out (default: `50000`).
    """)
    final Integer socketTimeout

    @ConfigOption
    @Description("""
        The target network throughput (in Gbps) used for S3 uploads and downloads (default: `10`).
    """)
    final Double targetThroughputInGbps

    // deprecated

    @Deprecated
    @ConfigOption
    @Description("""
        The size of a single part in a multipart upload (default: `100 MB`).
    """)
    final MemoryUnit uploadChunkSize

    @Deprecated
    @ConfigOption
    @Description("""
        The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).
    """)
    final Integer uploadMaxAttempts

    @Deprecated
    @ConfigOption
    @Description("""
        The maximum number of threads used for multipart upload (default: `10`).
    """)
    final Integer uploadMaxThreads

    @Deprecated
    @ConfigOption
    @Description("""
        The time to wait after a failed upload attempt to retry the part upload (default: `500ms`).
    """)
    final Duration uploadRetrySleep

    @Deprecated
    @ConfigOption
    @Description("""
        The S3 storage class applied to stored objects. Can be `STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, or `INTELLIGENT_TIERING` (default: `STANDARD`).
    """)
    final String uploadStorageClass

    private static final long _1MB = 1024 * 1024;
    // According to CRT Async client docs https://sdk.amazonaws.com/java/api/latest/software/amazon/awssdk/services/s3/S3CrtAsyncClientBuilder.html
    public static final long DEFAULT_PART_SIZE = 8 * _1MB;
    public static final int DEFAULT_INIT_BUFFER_PARTS = 10;
    // Maximum heap buffer size
    public static final long DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE = 400 * _1MB;

    AwsS3ClientConfig(Map opts) {
        super(opts)
        this.connectionTimeout = opts.connectionTimeout as Integer
        this.debug = opts.debug as Boolean
        if( endpoint && FileHelper.getUrlProtocol(endpoint) !in ['http','https'] )
            throw new IllegalArgumentException("S3 endpoint must begin with http:// or https:// prefix - offending value: '${endpoint}'")
        this.maxConcurrency = opts.maxConcurrency as Integer
        this.maxConnections = opts.maxConnections as Integer
        this.maxDownloadHeapMemory = opts.maxDownloadHeapMemory as MemoryUnit
        this.maxErrorRetry = opts.maxErrorRetry as Integer
        this.maxNativeMemory = opts.maxNativeMemory as MemoryUnit
        this.minimumPartSize = opts.minimumPartSize as MemoryUnit
        this.multipartThreshold = opts.multipartThreshold as MemoryUnit
        this.proxyHost = opts.proxyHost
        this.proxyPort = opts.proxyPort as Integer
        this.proxyScheme = opts.proxyScheme
        this.proxyUsername = opts.proxyUsername
        this.proxyPassword = opts.proxyPassword
        this.socketTimeout = opts.socketTimeout as Integer
        this.targetThroughputInGbps = opts.targetThroughputInGbps as Double
        this.uploadChunkSize = opts.uploadChunkSize as MemoryUnit
        this.uploadMaxAttempts = opts.uploadMaxAttempts as Integer
        this.uploadMaxThreads = opts.uploadMaxThreads as Integer
        this.uploadRetrySleep = opts.uploadRetrySleep as Duration
        checkDownloadBufferParams()
    }

    Map<String,String> getAwsClientConfig() {
        return [
            connection_timeout: connectionTimeout?.toString(),
            max_concurrency: maxConcurrency?.toString(),
            max_connections: maxConnections?.toString(),
            max_download_heap_memory: maxDownloadHeapMemory?.toBytes()?.toString(),
            max_error_retry: maxErrorRetry?.toString(),
            max_native_memory: maxNativeMemory?.toBytes()?.toString(),
            minimum_part_size: minimumPartSize?.toBytes()?.toString(),
            multipart_threshold: multipartThreshold?.toBytes()?.toString(),
            proxy_host: proxyHost?.toString(),
            proxy_port: proxyPort?.toString(),
            proxy_scheme: proxyScheme?.toString(),
            proxy_username: proxyUsername?.toString(),
            proxy_password: proxyPassword?.toString(),
            requester_pays: requesterPays?.toString(),
            s3_acl: s3Acl?.toString(),
            socket_timeout: socketTimeout?.toString(),
            storage_encryption: storageEncryption?.toString(),
            storage_kms_key_id: storageKmsKeyId?.toString(),
            target_throughput_in_gbps: targetThroughputInGbps?.toString(),
            upload_chunk_size: uploadChunkSize?.toBytes()?.toString(),
            upload_max_attempts: uploadMaxAttempts?.toString(),
            upload_max_threads: uploadMaxThreads?.toString(),
            upload_retry_sleep: uploadRetrySleep?.toMillis()?.toString(),
            upload_storage_class: storageClass?.toString()
        ].findAll { k, v -> v != null }
    }

    void checkDownloadBufferParams() {
        if( maxDownloadHeapMemory != null  && maxDownloadHeapMemory.toBytes() == 0L ) {
            throw new IllegalArgumentException("Configuration option `aws.client.maxDownloadHeapMemory` can't be 0")
        }
        if( minimumPartSize != null && minimumPartSize.toBytes() == 0L ) {
            throw new IllegalArgumentException("Configuration option `aws.client.minimumPartSize` can't be 0")
        }
        if( maxDownloadHeapMemory != null || minimumPartSize != null ) {
            final maxBuffer = maxDownloadHeapMemory ? maxDownloadHeapMemory.toBytes() : DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE
            final partSize = minimumPartSize ? minimumPartSize.toBytes() : DEFAULT_PART_SIZE
            if( maxBuffer < DEFAULT_INIT_BUFFER_PARTS * partSize ) {
                throw new IllegalArgumentException("Configuration option `aws.client.maxDownloadHeapMemory` must be at least " + DEFAULT_INIT_BUFFER_PARTS + " times `aws.client.minimumPartSize`")
            }
        }

    }
}
