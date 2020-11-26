/*
 * Copyright 2019, Google Inc
 * Copyright 2018, WuxiNextcode
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

package nextflow.cloud.google.lifesciences

import static nextflow.cloud.google.lifesciences.GoogleLifeSciencesHelper.*

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.CloudTransferOptions
import nextflow.executor.BashFunLib
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.Escape
/**
 * Defines the file/script copy strategies for Google Pipelines.
 *
 * @author Paolo Di Tommaso
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@Slf4j
@CompileStatic
class GoogleLifeSciencesFileCopyStrategy extends SimpleFileCopyStrategy {

    GoogleLifeSciencesConfig config
    GoogleLifeSciencesTaskHandler handler
    TaskBean task

    GoogleLifeSciencesFileCopyStrategy(TaskBean bean, GoogleLifeSciencesTaskHandler handler) {
        super(bean)
        this.handler = handler
        this.config = handler.executor.config
        this.task = bean
    }

    String getBeforeStartScript() {
        def maxConnect = config.maxParallelTransfers ?: CloudTransferOptions.MAX_TRANSFER
        def attempts = config.maxTransferAttempts ?: CloudTransferOptions.MAX_TRANSFER_ATTEMPTS
        def delayBetweenAttempts = config.delayBetweenAttempts ?: CloudTransferOptions.DEFAULT_DELAY_BETWEEN_ATTEMPTS

        BashFunLib.body(maxConnect, attempts, delayBetweenAttempts) +

        '''
        # google storage helper
        nxf_gs_download() {
            local source=$1
            local target=$2
            local project=$3
            local basedir=$(dirname $2)
            local ret
            local opts="-m -q"
            local opts
            if [[ $project ]]; then
              opts=('-q' '-m' '-u' "$project")
            else
              opts=('-q' '-m')
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
            gsutil -m -q cp -R "$name" "$target/$name"
        }
        '''.stripIndent()
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {

        final remoteTaskDir = getRemoteTaskDir(workDir)
        final stagingCommands = []

        final reqPay = config.enableRequesterPaysBuckets ? config.project : ''

        for( String target : inputFiles.keySet() ) {
            final sourcePath = inputFiles.get(target)

            final cmd = config.maxTransferAttempts > 1
                    ? "downloads+=(\"nxf_cp_retry nxf_gs_download ${Escape.uriPath(sourcePath)} ${Escape.path(target)} ${reqPay}\")"
                    : "downloads+=(\"nxf_gs_download ${Escape.uriPath(sourcePath)} ${Escape.path(target)} ${reqPay}\")"

            stagingCommands.add(cmd)
        }

        final result = new StringBuilder()
        // touch
        result.append("echo start | gsutil -q cp  -c - ${remoteTaskDir}/${TaskRun.CMD_START}").append('\n')

        // stage files
        if( stagingCommands ) {
            result.append('downloads=()\n')
            result.append(stagingCommands.join('\n')).append('\n')
            result.append('nxf_parallel "${downloads[@]}"\n')
        }

        // copy the remoteBinDir if it is defined
        if(config.remoteBinDir) {
            final localTaskDir = getLocalTaskDir(workDir)
            result
                    .append("mkdir -p $localTaskDir/nextflow-bin").append('\n')
                    .append("gsutil -m -q cp -P -r ${Escape.uriPath(config.remoteBinDir)}/* $localTaskDir/nextflow-bin").append('\n')
        }

        result.toString()
    }

    /**
     * Trim-trailing-slash
     * @return
     */
    private String tts(String loc) {
        while( log && loc.size()>1 && loc.endsWith('/') )
            loc = loc.substring(0,loc.length()-1)
        return loc
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        final patterns = normalizeGlobStarPaths(outputFiles)
        // create a bash script that will copy the out file to the working directory
        log.trace "[GLS] Unstaging file path: $patterns"

        if( !patterns )
            return null

        final escape = new ArrayList(outputFiles.size())
        for( String it : patterns )
            escape.add(tts(Escape.path(it)))

        """\
        uploads=()
        IFS=\$'\\n'
        for name in \$(eval "ls -1d ${escape.join(' ')}" | sort | uniq); do
            uploads+=("nxf_gs_upload '\$name' ${Escape.uriPath(targetDir)}")
        done
        unset IFS
        nxf_parallel "\${uploads[@]}"
        """.stripIndent()
    }


    String copyFile(String local, Path target) {
        /*
         * -m = run in parallel
         * -q = quiet mode
         * cp = copy
         * -R = recursive copy
         */
        "gsutil -m -q cp -R ${Escape.path(local)} ${Escape.uriPath(target)}"
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, boolean container) {
        if( container )
            throw new IllegalArgumentException("Parameter `container` not supported by ${this.class.simpleName}")

        final localTaskDir = getLocalTaskDir(workDir)
        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        // remove any external PATH
        if( copy.containsKey('PATH') )
            copy.remove('PATH')
        // when a remote bin directory is provide managed it properly
        if( config.remoteBinDir ) {
            result << "chmod +x $localTaskDir/nextflow-bin/*\n"
            result << "export PATH=$localTaskDir/nextflow-bin:\$PATH\n"
        }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,false)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }

}
