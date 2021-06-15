/*
 * Copyright 2020-2021, Seqera Labs
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

package nextflow.cloud.azure.file

import groovy.transform.Memoized
import nextflow.executor.BashFunLib
import nextflow.util.Duration

/**
 * Azure Bash helper functions
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBashLib extends BashFunLib<AzBashLib> {

    protected String azLib() {
        '''
        nxf_az_upload() {
            local name=$1
            local target=${2%/} ## remove ending slash
        
            if [[ -d $name ]]; then
              azcopy cp "$name" "$target?$AZ_SAS" --recursive
            else
              azcopy cp "$name" "$target/$name?$AZ_SAS"
            fi
        }
        
        nxf_az_download() {
            local source=$1
            local target=$2
            local basedir=$(dirname $2)
            local ret
            mkdir -p "$basedir"
        
            ret=$(azcopy cp "$source?$AZ_SAS" "$target" 2>&1) || {
                ## if fails check if it was trying to download a directory
                mkdir -p $target
                azcopy cp "$source/*?$AZ_SAS" "$target" --recursive >/dev/null || {
                    rm -rf $target
                    >&2 echo "Unable to download path: $source"
                    exit 1
                }
            }
        }
        '''.stripIndent()
    }

    String render() {
        super.render() + azLib()
    }

    @Memoized
    static String script(Integer maxParallelTransfers, Integer maxTransferAttempts, Duration delayBetweenAttempts) {
        new AzBashLib()
                .includeCoreFun(true)
                .withMaxParallelTransfers(maxParallelTransfers)
                .withMaxTransferAttempts(maxTransferAttempts)
                .withDelayBetweenAttempts(delayBetweenAttempts)
                .render()
    }

    @Memoized
    static String script() {
        new AzBashLib().render()
    }
}
