/*
 * Copyright 2020-2025, Seqera Labs
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
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.file.FileHelper
import nextflow.script.dsl.Description
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.s3.model.ObjectCannedACL

import static nextflow.cloud.aws.util.AwsHelper.parseS3Acl

/**
 * Model AWS S3 bucket config settings
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class AwsS3CommonConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        Allow the access of public S3 buckets without providing AWS credentials (default: `false`). Any service that does not accept unsigned requests will return a service access error.
    """)
    final Boolean anonymous

    @ConfigOption
    @Description("""
        The AWS S3 API entry point e.g. `https://s3-us-west-1.amazonaws.com`. The endpoint must include the protocol prefix e.g. `https://`.
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        S3 Bucket specific AWS region (e.g. `us-east-1`).
    """)
    final String region

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

    AwsS3CommonConfig(Map opts) {
        this.anonymous = opts.anonymous as Boolean
        this.endpoint = opts.endpoint ?: SysEnv.get('AWS_S3_ENDPOINT')
        if( endpoint && FileHelper.getUrlProtocol(endpoint) !in ['http','https'] )
            throw new IllegalArgumentException("S3 endpoint must begin with http:// or https:// prefix - offending value: '${endpoint}'")
        this.region = opts.region as String
        this.requesterPays = opts.requesterPays as Boolean
        this.s3Acl = parseS3Acl(opts.s3Acl as String)
        this.s3PathStyleAccess = opts.s3PathStyleAccess as Boolean
        this.storageClass = parseStorageClass((opts.storageClass ?: opts.uploadStorageClass) as String)     // 'uploadStorageClass' is kept for legacy purposes
        this.storageEncryption = parseStorageEncryption(opts.storageEncryption as String)
        this.storageKmsKeyId = opts.storageKmsKeyId as String
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

    /**
     * Looks for the region defined in endpoints such as https://xxx.<region>.amazonaws.com
     * @returns Region defined in the endpoint. Null if no endpoint or custom endpoint is defined,
     * or when URI region subdomain doesn't match with a region (global or multi-region access point)
     */
    String getEndpointRegion(){
        if( !endpoint || isCustomEndpoint() )
            return null

        try {
            String host = URI.create(endpoint).getHost()
            final hostDomains = host.split('\\.')
            if (hostDomains.size() < 3) {
                log.debug("Region subdomain doesn't exist in endpoint '${endpoint}'")
                return null
            }
            final region = hostDomains[hostDomains.size()-3]
            if (!Region.regions().contains(Region.of(region))){
                log.debug("Region '${region}' extracted from endpoint '${endpoint}' is not valid")
                return null
            }
            return region

        } catch (Exception e){
            log.debug("Exception getting region from endpoint: '${endpoint}' - ${e.message}")
            return null
        }
    }

    Map<String, Object> toBucketConfigMap(){
        return ([
            anonymous: anonymous,
            endpoint: endpoint,
            requesterPays: requesterPays,
            region: region,
            s3PathStyleAccess: s3PathStyleAccess,
            s3Acl: s3Acl?.toString(),
            storageEncryption: storageEncryption,
            storageKmsKeyId: storageKmsKeyId,
            storageClass: storageClass
        ] as Map<String, Object>).findAll { k, v -> v != null }
    }

    Map<String,String> toLegacyConfig() {
        return [
            requester_pays: requesterPays?.toString(),
            s3_acl: s3Acl?.toString(),
            storage_encryption: storageEncryption?.toString(),
            storage_kms_key_id: storageKmsKeyId?.toString(),
            upload_storage_class: storageClass?.toString()
        ].findAll { k, v -> v != null }
    }

}
