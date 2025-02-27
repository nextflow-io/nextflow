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
 */

package nextflow.cloud.google.util

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Session
import nextflow.cloud.google.lifesciences.GoogleLifeSciencesConfig
import nextflow.executor.BashFunLib

/**
 * Provide Bash helpers for Google Storage handling
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GsBashLib extends BashFunLib<GsBashLib> {

    private int parallelThreadCount = GoogleLifeSciencesConfig.DEF_PARALLEL_THREAD_COUNT
    private int downloadMaxComponents = GoogleLifeSciencesConfig.DEF_DOWNLOAD_MAX_COMPONENTS

    GsBashLib withParallelThreadCount(Integer value) {
        if( value )
            parallelThreadCount = value
        return this
    }

    GsBashLib withDownloadMaxComponents(Integer value) {
        if( value )
            downloadMaxComponents = value
        return this
    }

    protected makeEnv() {
        """
        # google storage helper
        gs_opts=('-q' '-m' '-o' 'GSUtil:parallel_thread_count=$parallelThreadCount' '-o' 'GSUtil:sliced_object_download_max_components=$downloadMaxComponents')
        """.stripIndent(true)
    }

    protected makeLib() {
        '''
        nxf_gs_download() {
            local source=$1
            local target=$2
            local project=$3
            local basedir=$(dirname $2)
            local ret
            local opts=("${gs_opts[@]}")
            if [[ $project ]]; then
              opts+=('-u' "$project")
            fi
             
            ## download assuming it's a file download
            mkdir -p $basedir
            ret=$(gsutil ${opts[@]} cp "$source" "$target" 2>&1) || {
                ## if fails check if it was trying to download a directory
                mkdir $target
                gsutil ${opts[@]} cp -R "$source/*" "$target" || {
                  rm -rf $target
                  >&2 echo "Unable to download path: $source"
                  exit 1
                }
            }
        }

        nxf_gs_upload() {
            local name=$1
            local target=$2
            gsutil ${gs_opts[@]} cp -R "$name" "$target/$name"
        }
        '''.stripIndent(true)
    }

    @Override
    String render() {
        super.render() + makeEnv() + makeLib()
    }

    @Memoized
    static String fromConfig(GoogleLifeSciencesConfig config) {
        new GsBashLib()
                .includeCoreFun(true)
                .withMaxParallelTransfers(config.maxParallelTransfers)
                .withMaxTransferAttempts(config.maxTransferAttempts)
                .withDelayBetweenAttempts(config.delayBetweenAttempts)
                .withParallelThreadCount(config.parallelThreadCount)
                .withDownloadMaxComponents(config.downloadMaxComponents)
                .render()
    }

    @Memoized
    static String fromSession(Session session) {
        final parallelThreadCount = session.config.navigate('google.storage.parallelThreadCount') as Integer
        final downloadMaxComponents = session.config.navigate('google.storage.downloadMaxComponents') as Integer

        new GsBashLib()
                .withParallelThreadCount(parallelThreadCount)
                .withDownloadMaxComponents(downloadMaxComponents)
                .render()
    }
}
