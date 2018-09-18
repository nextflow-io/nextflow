/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
        if( value in [null, 'STANDARD', 'REDUCED_REDUNDANCY'])
            this.storageClass = value
        else
            log.warn "Unsupported AWS storage-class: $value"
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