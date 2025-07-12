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
        Allow the access of public S3 buckets without providing AWS credentials. Any service that does not accept unsigned requests will return a service access error.
    """)
    final Boolean anonymous

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) when initially establishing a connection before timing out.
    """)
    final int connectionTimeout

    final Boolean debug

    @ConfigOption
    @Description("""
        The AWS S3 API entry point e.g. `https://s3-us-west-1.amazonaws.com`. The endpoint must include the protocol prefix e.g. `https://`.
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        The maximum number of concurrency in S3 async clients.
    """)
    final Integer maxConcurrency

    @ConfigOption
    @Description("""
        The maximum number of allowed open HTTP connections.
    """)
    final Integer maxConnections

    @ConfigOption
    @Description("""
        The maximum number of retry attempts for failed retryable requests.
    """)
    final int maxErrorRetry

    @ConfigOption
    @Description("""
        The maximum native memory used by the S3 asynchronous client for S3 transfers.
    """)
    final MemoryUnit maxNativeMemory

    @ConfigOption
    @Description("""
        The minimum size of a single part in a multipart upload (default: `8 MB`).
    """)
    final MemoryUnit minimumPartSize

    @ConfigOption
    @Description("""
        The S3 Async client threshold to create multipart S3 transfers. Default is the same as `minimumPartSize`.
    """)
    final MemoryUnit multipartThreshold

    @ConfigOption
    @Description("""
        The proxy host to connect through.
    """)
    final String proxyHost

    @ConfigOption
    @Description("""
        The port on the proxy host to connect through.
    """)
    final int proxyPort

    @ConfigOption
    @Description("""
        The protocol scheme to use when connecting through a proxy (http/https).
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
        Enable the requester pays feature for S3 buckets.
    """)
    final Boolean requesterPays

    @ConfigOption(types=[String])
    @Description("""
        Specify predefined bucket permissions, also known as *canned ACL*. Can be one of `Private`, `PublicRead`, `PublicReadWrite`, `AuthenticatedRead`, `LogDeliveryWrite`, `BucketOwnerRead`, `BucketOwnerFullControl`, or `AwsExecRead`.
        
        [Read more](https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl)
    """)
    final ObjectCannedACL s3Acl

    @ConfigOption
    @Description("""
        Enable the use of path-based access model that is used to specify the address of an object in S3-compatible storage systems.
    """)
    final Boolean s3PathStyleAccess

    @ConfigOption
    @Description("""
        The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out.
    """)
    final int socketTimeout

    @ConfigOption
    @Description("""
        The S3 storage class applied to stored objects, one of \\[`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`\\] (default: `STANDARD`).
    """)
    final String storageClass

    @ConfigOption
    @Description("""
        The S3 server side encryption to be used when saving objects on S3, either `AES256` or `aws:kms` values are allowed.
    """)
    final String storageEncryption

    @ConfigOption
    @Description("""
        The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket.
    """)
    final String storageKmsKeyId

    @ConfigOption
    @Description("""
        The S3 Async client target network throughput in Gbps. This value is used to automatically set `maxConcurrency` and `maxNativeMemory` (default: `10`).
    """)
    final Double targetThroughputInGbps

    @ConfigOption
    @Description("""
        The number of threads used by the S3 transfer manager (default: `10`).
    """)
    final Integer transferManagerThreads

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
        The maximum number of threads used for multipart upload.
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
        The S3 storage class applied to stored objects, one of \\[`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`\\] (default: `STANDARD`).
    """)
    final String uploadStorageClass

    AwsS3Config(Map opts) {
        this.anonymous = opts.anonymous as Boolean
        this.debug = opts.debug as Boolean
        this.endpoint = opts.endpoint ?: SysEnv.get('AWS_S3_ENDPOINT')
        if( endpoint && FileHelper.getUrlProtocol(endpoint) !in ['http','https'] )
            throw new IllegalArgumentException("S3 endpoint must begin with http:// or https:// prefix - offending value: '${endpoint}'")
        this.maxConcurrency = opts.maxConcurrency as Integer
        this.maxConnections = opts.maxConnections as Integer
        this.maxNativeMemory = opts.maxNativeMemory as MemoryUnit
        this.minimumPartSize = opts.minimumPartSize as MemoryUnit
        this.multipartThreshold = opts.multipartThreshold as MemoryUnit
        this.requesterPays = opts.requesterPays as Boolean
        this.s3Acl = parseS3Acl(opts.s3Acl as String)
        this.s3PathStyleAccess = opts.s3PathStyleAccess as Boolean
        this.storageClass = parseStorageClass((opts.storageClass ?: opts.uploadStorageClass) as String)     // 'uploadStorageClass' is kept for legacy purposes
        this.storageEncryption = parseStorageEncryption(opts.storageEncryption as String)
        this.storageKmsKeyId = opts.storageKmsKeyId
        this.targetThroughputInGbps = opts.targetThroughputInGbps as Double
        this.transferManagerThreads = opts.transferManagerThreads as Integer
        this.uploadChunkSize = opts.uploadChunkSize as MemoryUnit
        this.uploadMaxAttempts = opts.uploadMaxAttempts as Integer
        this.uploadMaxThreads = opts.uploadMaxThreads as Integer
        this.uploadRetrySleep = opts.uploadRetrySleep as Duration
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
        return endpoint && !endpoint.endsWith(".amazonaws.com")
    }

    Map<String,?> getAwsClientConfig() {
        return [
            max_concurrency: maxConcurrency,
            max_native_memory: maxNativeMemory,
            minimum_part_size: minimumPartSize?.toBytes(),
            multipart_threshold: multipartThreshold?.toBytes(),
            target_throughput_in_gbps: targetThroughputInGbps,
            transfer_manager_threads: transferManagerThreads,
            upload_chunk_size: uploadChunkSize?.toBytes(),
            upload_max_attempts: uploadMaxAttempts,
            upload_max_threads: uploadMaxThreads,
            upload_retry_sleep: uploadRetrySleep?.toMillis(),
        ]
    }
}
