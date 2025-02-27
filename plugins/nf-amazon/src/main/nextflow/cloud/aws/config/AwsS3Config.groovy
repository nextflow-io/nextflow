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

import com.amazonaws.services.s3.model.CannedAccessControlList
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.file.FileHelper
/**
 * Model AWS S3 config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsS3Config {

    private String endpoint

    private String storageClass

    private String storageEncryption

    private String storageKmsKeyId

    private Boolean debug

    private CannedAccessControlList s3Acl

    private Boolean pathStyleAccess

    private Boolean anonymous

    private Boolean requesterPays

    AwsS3Config(Map opts) {
        this.debug = opts.debug as Boolean
        this.endpoint = opts.endpoint ?: SysEnv.get('AWS_S3_ENDPOINT')
        if( endpoint && FileHelper.getUrlProtocol(endpoint) !in ['http','https'] )
            throw new IllegalArgumentException("S3 endpoint must begin with http:// or https:// prefix - offending value: '${endpoint}'")
        this.storageClass = parseStorageClass((opts.storageClass ?: opts.uploadStorageClass) as String)     // 'uploadStorageClass' is kept for legacy purposes
        this.storageEncryption = parseStorageEncryption(opts.storageEncryption as String)
        this.storageKmsKeyId = opts.storageKmsKeyId
        this.pathStyleAccess = opts.s3PathStyleAccess as Boolean
        this.anonymous = opts.anonymous as Boolean
        this.s3Acl = parseS3Acl(opts.s3Acl as String)
        this.requesterPays = opts.requesterPays as Boolean
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
    String getEndpoint() {
        return endpoint
    }

    String getStorageClass() {
        return storageClass
    }

    String getStorageEncryption() {
        return storageEncryption
    }

    String getStorageKmsKeyId() {
        return storageKmsKeyId
    }

    Boolean getDebug() {
        return debug
    }

    CannedAccessControlList getS3Acl() {
        return s3Acl
    }

    Boolean getPathStyleAccess() {
        return pathStyleAccess
    }

    Boolean getAnonymous() {
        return anonymous
    }

    Boolean getRequesterPays() {
        return requesterPays
    }

    boolean isCustomEndpoint() {
        endpoint && !endpoint.contains(".amazonaws.com")
    }
}
