/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.executor

import groovy.transform.CompileStatic
import nextflow.cloud.CloudTransferOptions
import nextflow.util.Duration

/**
 * Bash common functions library
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BashFunLib<V extends BashFunLib> {

    protected int maxParallelTransfers = CloudTransferOptions.MAX_TRANSFER
    protected int maxTransferAttempts = CloudTransferOptions.MAX_TRANSFER_ATTEMPTS
    protected Duration delayBetweenAttempts = CloudTransferOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS
    protected boolean includeCoreFun

    V includeCoreFun(boolean value) {
        this.includeCoreFun = value
        return (V)this
    }

    V withMaxParallelTransfers(Integer value) {
        if( value )
            maxParallelTransfers = value
        return (V)this
    }

    V withMaxTransferAttempts(Integer value) {
        if( value )
            maxTransferAttempts = value
        return (V)this
    }

    V withDelayBetweenAttempts(Duration value) {
        if( value )
            delayBetweenAttempts = value
        return (V)this
    }

    String coreLib() {
        """\
        # bash helper functions
        nxf_cp_retry() {
            local max_attempts=$maxTransferAttempts
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
              nxf_sleep \$timeout
              attempt=\$(( attempt + 1 ))
              timeout=\$(( timeout * 2 ))
            done
        }
        
        nxf_parallel() {
            IFS=\$'\\n'
            local cmd=("\$@")
            local cpus=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
            local max=\$(if (( cpus>$maxParallelTransfers )); then echo $maxParallelTransfers; else echo \$cpus; fi)
            local i=0
            local pid=()
            (
            set +u
            while ((i<\${#cmd[@]})); do
                local copy=()
                for x in "\${pid[@]}"; do
                  # if the process exist, keep in the 'copy' array, otherwise wait on it to capture the exit code
                  # see https://github.com/nextflow-io/nextflow/pull/4050
                  [[ -e /proc/\$x ]] && copy+=(\$x) || wait \$x
                done
                pid=("\${copy[@]}")
                
                if ((\${#pid[@]}>=\$max)); then
                  nxf_sleep 0.2
                else
                  eval "\${cmd[\$i]}" &
                  pid+=(\$!)
                  ((i+=1))
                fi
            done
            for p in "\${pid[@]}"; do
                wait \$p
            done
            )
            unset IFS
        }
        """.stripIndent(true)
    }

    String render() {
        includeCoreFun ? coreLib() : ''
    }

    static String body(int maxParallelTransfer, int maxTransferAttempts, Duration delay) {
        new BashFunLib()
                .includeCoreFun(true)
                .withMaxParallelTransfers(maxParallelTransfer)
                .withMaxTransferAttempts(maxTransferAttempts)
                .withDelayBetweenAttempts(delay)
                .render()
    }
}
