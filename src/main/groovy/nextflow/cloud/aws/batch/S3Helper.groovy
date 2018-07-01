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

/**
 * AWS S3 helper class
 */
class S3Helper {

    static String getUploaderScript(AwsOptions opts) {
        def cli = opts.getAwsCli()
        def storage = opts.storageClass ?: 'STANDARD'
        def encryption = opts.storageEncryption ? "--sse $opts.storageEncryption " : ''

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
            local cmd=("\$@")
            local cpus=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
            local max=\$(if (( cpus>16 )); then echo 16; else echo \$cpus; fi)
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
        }     
        """.stripIndent()
    }

}