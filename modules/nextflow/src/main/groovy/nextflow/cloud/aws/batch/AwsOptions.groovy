/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
 */

package nextflow.cloud.aws.batch

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException

/**
 * Helper class wrapping AWS config options required for Batch job executions
 */
@Slf4j
@ToString
@EqualsAndHashCode
@CompileStatic
class AwsOptions {

    String cliPath

    String storageClass

    String storageEncryption

    String remoteBinDir

    String region

    AwsOptions() { }

    AwsOptions(Session session) {
        storageClass = session.config.navigate('aws.client.uploadStorageClass') as String
        storageEncryption = session.config.navigate('aws.client.storageEncryption') as String
        region = session.config.navigate('aws.region') as String
    }

    void setStorageClass(String value) {
        if( value in [null, 'STANDARD', 'STANDARD_IA', 'ONEZONE_IA', 'INTELLIGENT_TIERING', 'REDUCED_REDUNDANCY' ]) {
            this.storageClass = value
            if (value == 'REDUCED_REDUNDANCY') {
                log.warn "AWS S3 Storage Class `REDUCED_REDUNDANCY` is deprecated (and more expensive than `STANDARD`). For cost savings, look to `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`."
            }
        } else {
            log.warn "Unsupported AWS storage-class: $value"
        }
    }

    void setStorageEncryption(String value) {
        if( value in [null,'AES256'] )
            this.storageEncryption = value
        else
            log.warn "Unsupported AWS storage-encryption: $value"
    }

    void setCliPath(String value) {
        if( !value )
            this.cliPath = null
        else {
            if( !value.startsWith('/') ) throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must be an absolute path")
            if( !value.endsWith('/bin/aws')) throw new ProcessUnrecoverableException("Not a valid aws-cli tools path: $value -- it must end with the `/bin/aws` suffix")
            this.cliPath = value
        }
    }

    String getAwsCli() {
        def result = getCliPath()
        if( !result ) result = 'aws'
        if( region ) result += " --region $region"
        return result
    }
}