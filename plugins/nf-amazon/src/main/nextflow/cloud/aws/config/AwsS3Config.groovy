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

import static nextflow.cloud.aws.util.AwsHelper.*

import software.amazon.awssdk.services.s3.model.ObjectCannedACL
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.script.dsl.Description
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
/**
 * Model AWS S3 config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsS3Config implements ConfigScope {

    @ConfigOption
    @Description("""
        Allow the access of public S3 buckets without providing AWS credentials (default: `false`). Any service that does not accept unsigned requests will return a service access error.
    """)
    final Boolean anonymous

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) when initially establishing a connection before timing out (default: `10000`).
    """)
    final Integer connectionTimeout

    final Boolean debug

    @ConfigOption
    @Description("""
        The AWS S3 API entry point e.g. `https://s3-us-west-1.amazonaws.com`. The endpoint must include the protocol prefix e.g. `https://`.
    """)
    final String endpoint

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
        Use [Requester Pays](https://docs.aws.amazon.com/AmazonS3/latest/userguide/RequesterPaysBuckets.html) for S3 buckets (default: `false`).
    """)
    final Boolean requesterPays

    @ConfigOption(types=[String])
    @Description("""
        Specify predefined bucket permissions, also known as [canned ACL](https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl). Can be one of `Private`, `PublicRead`, `PublicReadWrite`, `AuthenticatedRead`, `LogDeliveryWrite`, `BucketOwnerRead`, `BucketOwnerFullControl`, or `AwsExecRead`.
    """)
    final ObjectCannedACL s3Acl

    @ConfigOption
    @Description("""
        Use the path-based access model to access objects in S3-compatible storage systems (default: `false`).
    """)
    final Boolean s3PathStyleAccess

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out (default: `50000`).
    """)
    final Integer socketTimeout

    @ConfigOption
    @Description("""
        The S3 storage class applied to stored objects, one of \\[`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`\\] (default: `STANDARD`).
    """)
    final String storageClass

    @ConfigOption
    @Description("""
        The S3 server side encryption to be used when saving objects on S3. Can be `AES256` or `aws:kms` (default: none).
    """)
    final String storageEncryption

    @ConfigOption
    @Description("""
        The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket.
    """)
    final String storageKmsKeyId

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

    // According to CRT Async client docs https://sdk.amazonaws.com/java/api/latest/software/amazon/awssdk/services/s3/S3CrtAsyncClientBuilder.html
    public static final long DEFAULT_PART_SIZE = MemoryUnit.of('8 MB').toBytes()
    public static final int DEFAULT_INIT_BUFFER_PARTS = 10
    // Maximum heap buffer size
    public static final long DEFAULT_MAX_DOWNLOAD_BUFFER_SIZE = MemoryUnit.of('400 MB').toBytes()

    AwsS3Config(Map opts) {
        this.anonymous = opts.anonymous as Boolean
        this.connectionTimeout = opts.connectionTimeout as Integer
        this.debug = opts.debug as Boolean
        this.endpoint = opts.endpoint ?: SysEnv.get('AWS_S3_ENDPOINT')
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
        this.requesterPays = opts.requesterPays as Boolean
        this.s3Acl = parseS3Acl(opts.s3Acl as String)
        this.s3PathStyleAccess = opts.s3PathStyleAccess as Boolean
        this.socketTimeout = opts.socketTimeout as Integer
        this.storageClass = parseStorageClass((opts.storageClass ?: opts.uploadStorageClass) as String)     // 'uploadStorageClass' is kept for legacy purposes
        this.storageEncryption = parseStorageEncryption(opts.storageEncryption as String)
        this.storageKmsKeyId = opts.storageKmsKeyId
        this.targetThroughputInGbps = opts.targetThroughputInGbps as Double
        this.uploadChunkSize = opts.uploadChunkSize as MemoryUnit
        this.uploadMaxAttempts = opts.uploadMaxAttempts as Integer
        this.uploadMaxThreads = opts.uploadMaxThreads as Integer
        this.uploadRetrySleep = opts.uploadRetrySleep as Duration
        checkDownloadBufferParams()
    }

    private String parseStorageClass(String value) {
        if( value in [null, 'STANDARD', 'STANDARD_IA', 'ONEZONE_IA', 'INTELLIGENT_TIERING', 'REDUCED_REDUNDANCY' ]) {
            if (value == 'REDUCED_REDUNDANCY') {
                log.warn "AWS S3 Storage Class `REDUCED_REDUNDANCY` is deprecated (and more expensive than `STANDARD`). For cost savings, look to `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`."
            }
            return value
        } else {
            log.warn "Unsupported AWS storage-class: $value"
            return null
        }
    }

    private String parseStorageEncryption(String value) {
        if( value in [null,'AES256','aws:kms'] )
            return value
        //
        log.warn "Unsupported AWS storage-encryption: $value"
        return null
    }

    // ==== getters =====

    Boolean getPathStyleAccess() {
        return s3PathStyleAccess
    }

    boolean isCustomEndpoint() {
        endpoint && !endpoint.endsWith(".amazonaws.com")
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
