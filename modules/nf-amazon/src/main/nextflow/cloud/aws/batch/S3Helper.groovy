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

import nextflow.Global
import nextflow.Session

/**
 * AWS S3 helper class
 */
class S3Helper {

    static String getUploaderScript(AwsOptions opts) {
        def cli = opts.getAwsCli()
        def storage = opts.storageClass ?: 'STANDARD'
        def encryption = opts.storageEncryption ? "--sse $opts.storageEncryption " : ''
        def maxConnect = opts.maxParallelTransfers ?: AwsOptions.MAX_TRANSFER
        def attempts = opts.maxTransferAttempts ?: AwsOptions.MAX_TRANSFER_ATTEMPTS
        def delayBetweenAttempts = opts.delayBetweenAttempts ?: AwsOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS

        """
        # aws helper
        nxf_s3_upload() {
            local pattern=\$1
            local s3path=\$2
            IFS=\$'\\n'
            for name in \$(eval "ls -1d \$pattern");do
              if [[ -d "\$name" ]]; then
                $cli s3 cp --only-show-errors --recursive $encryption--storage-class $storage "\$name" "\$s3path/\$name"
              else
                $cli s3 cp --only-show-errors $encryption--storage-class $storage "\$name" "\$s3path/\$name"
              fi
            done
            unset IFS
        }
        
        nxf_s3_retry() {
            local max_attempts=$attempts
            local timeout=${delayBetweenAttempts.seconds}
            local attempt=0
            local exitCode=0
            while (( \$attempt < \$max_attempts ))
            do
              if "\$@"
                then
                  return 0
              else
                exitCode=\$?
              fi
              if [[ \$exitCode == 0 ]]
              then
                break
              fi
              sleep \$timeout
              attempt=\$(( attempt + 1 ))
              timeout=\$(( timeout * 2 ))
            done
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
     
        nxf_parallel() {
            IFS=\$'\\n'
            local cmd=("\$@")
            local cpus=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
            local max=\$(if (( cpus>$maxConnect )); then echo $maxConnect; else echo \$cpus; fi)
            local i=0
            local pid=()
            (
            set +u
            while ((i<\${#cmd[@]})); do
                local copy=()
                for x in "\${pid[@]}"; do
                  [[ -e /proc/\$x ]] && copy+=(\$x) 
                done
                pid=("\${copy[@]}")
                
                if ((\${#pid[@]}>=\$max)); then 
                  sleep 1 
                else 
                  eval "\${cmd[\$i]}" &
                  pid+=(\$!)
                  ((i+=1))
                fi 
            done
            ((\${#pid[@]}>0)) && wait \${pid[@]}
            )
            unset IFS
        }
        """.stripIndent()
    }

    static String getUploaderScript() {
        final opts = new AwsOptions(Global.session as Session)
        getUploaderScript(opts)
    }

}
