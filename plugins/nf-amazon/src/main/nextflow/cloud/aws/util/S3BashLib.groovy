/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cloud.aws.util

import nextflow.Global
import nextflow.Session
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.executor.BashFunLib

/**
 * AWS S3 helper class
 */
class S3BashLib extends BashFunLib<S3BashLib> {

    private String storageClass = 'STANDARD'
    private String storageEncryption = ''
    private String storageKmsKeyId = ''
    private String debug = ''
    private String cli = 'aws'
    private String retryMode

    S3BashLib withCliPath(String cliPath) {
        if( cliPath )
            this.cli = cliPath
        return this
    }

    S3BashLib withRetryMode(String value) {
        if( value )
            retryMode = value
        return this
    }

    S3BashLib withDebug(Boolean  value) {
        this.debug = value ? '--debug ' : ''
        return this
    }

    S3BashLib withStorageClass(String value) {
        if( value )
            this.storageClass = value
        return this
    }

    S3BashLib withStorageEncryption(String value) {
        if( value )
            this.storageEncryption = value ? "--sse $value " : ''
        return this
    }

    S3BashLib withStorageKmsKeyId(String value) {
        if( value )
            this.storageKmsKeyId = value ? "--sse-kms-key-id $value " : ''
        return this
    }

    protected String retryEnv() {
        if( !retryMode )
            return ''
        """
        # aws cli retry config
        export AWS_RETRY_MODE=${retryMode} 
        export AWS_MAX_ATTEMPTS=${maxTransferAttempts}
        """.stripIndent().rightTrim()
    }

    protected String s3Lib() {
        """
        # aws helper
        nxf_s3_upload() {
            local name=\$1
            local s3path=\$2
            if [[ "\$name" == - ]]; then
              $cli s3 cp --only-show-errors ${debug}${storageEncryption}${storageKmsKeyId}--storage-class $storageClass - "\$s3path"
            elif [[ -d "\$name" ]]; then
              $cli s3 cp --only-show-errors --recursive ${debug}${storageEncryption}${storageKmsKeyId}--storage-class $storageClass "\$name" "\$s3path/\$name"
            else
              $cli s3 cp --only-show-errors ${debug}${storageEncryption}${storageKmsKeyId}--storage-class $storageClass "\$name" "\$s3path/\$name"
            fi
        }
        
        nxf_s3_download() {
            local source=\$1
            local target=\$2
            local file_name=\$(basename \$1)
            local is_dir=\$($cli s3 ls \$source | grep -F "PRE \${file_name}/" -c)
            if [[ \$is_dir == 1 ]]; then
                $cli s3 cp --only-show-errors --recursive "\$source" "\$target"
            else 
                $cli s3 cp --only-show-errors "\$source" "\$target"
            fi
        }
        """.stripIndent()
    }

    String render() {
        super.render() + retryEnv() + s3Lib()
    }

    static private S3BashLib lib0(AwsOptions opts, boolean includeCore) {
        new S3BashLib()
                .includeCoreFun(includeCore)
                .withMaxParallelTransfers( opts.maxParallelTransfers )
                .withDelayBetweenAttempts(opts.delayBetweenAttempts )
                .withMaxTransferAttempts( opts.maxTransferAttempts )
                .withCliPath( opts.awsCli )
                .withStorageClass(opts.storageClass )
                .withStorageEncryption( opts.storageEncryption )
                .withStorageKmsKeyId( opts.storageKmsKeyId )
                .withRetryMode( opts.retryMode )
                .withDebug( opts.debug )
    }

    static String script(AwsOptions opts) {
        lib0(opts,true).render()
    }

    static String script() {
        final opts = new AwsOptions(Global.session as Session)
        lib0(opts,false).render()
    }
}
