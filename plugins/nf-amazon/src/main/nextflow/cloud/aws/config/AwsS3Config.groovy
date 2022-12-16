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
import nextflow.SysEnv
/**
 * Model AWS S3 config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AwsS3Config {

    private String endpoint

    private String storageClass

    private String storageEncryption

    private String storageKmsKeyId

    private Boolean debug

    AwsS3Config(Map opts) {
        this.debug = opts.debug as Boolean
        storageClass = opts.storageClass ?: opts.uploadStorageClass     // 'uploadStorageClass' is kept for legacy purposes
        this.endpoint = opts.endpoint ?: SysEnv.get('AWS_S3_ENDPOINT')
        this.storageKmsKeyId = opts.storageKmsKeyId
        this.storageEncryption = opts.storageEncryption
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

}
