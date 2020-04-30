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

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.Escape
import static nextflow.cloud.google.lifesciences.GoogleLifeSciencesHelper.getLocalTaskDir
import static nextflow.cloud.google.lifesciences.GoogleLifeSciencesHelper.getRemoteTaskDir


/**
 * Defines the file/script copy strategies for Google Pipelines.
 *
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

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {

        final localTaskDir = getLocalTaskDir(workDir)
        final remoteTaskDir = getRemoteTaskDir(workDir)
        final createDirectories  = []
        final stagingCommands = []

        final gsutilPrefix = new StringBuilder()
        gsutilPrefix.append("gsutil -m -q")

        if(config.enableRequesterPaysBuckets) {
            gsutilPrefix.append(" -u ${config.project}")
        }

        for( String stageName : inputFiles.keySet() ) {
            final storePath = inputFiles.get(stageName)
            final storePathIsDir = storePath.isDirectory()
            final stagePath = Paths.get(stageName)
            final parent = stagePath.parent
            final escapedStoreUri = Escape.uriPath(storePath)
            final escapedStageName = Escape.path(stageName)

            //check if we need to create parent dirs for staging file since gsutil doesn't create them for us
            if(parent) {
                createDirectories << "mkdir -p $localTaskDir/${Escape.path(parent)}".toString()
            }

            if(storePathIsDir) {
                stagingCommands << "$gsutilPrefix cp -R $escapedStoreUri/ $localTaskDir".toString()
                //check if we need to move the directory (gsutil doesn't support renaming directories on copy)
                if(parent || !storePath.toString().endsWith(stageName)) {
                    stagingCommands << "mv $localTaskDir/${Escape.path(storePath.name)} $localTaskDir/$escapedStageName".toString()
                }
            } else {
                stagingCommands << "$gsutilPrefix cp $escapedStoreUri $localTaskDir/$escapedStageName".toString()
            }
        }

        final result = new StringBuilder()
        // touch
        result.append("echo start | gsutil -q cp  -c - ${remoteTaskDir}/${TaskRun.CMD_START}").append('\n')

        // create directories
        if( createDirectories )
            result.append(createDirectories.join('\n')) .append('\n')

        // stage files
        if( stagingCommands )
            result.append(stagingCommands.join('\n')) .append('\n')

        // copy the remoteBinDir if it is defined
        if(config.remoteBinDir) {
            result
                    .append("mkdir -p $localTaskDir/nextflow-bin").append('\n')
                    .append("gsutil -m -q cp -P -r ${Escape.uriPath(config.remoteBinDir)}/* $localTaskDir/nextflow-bin").append('\n')
        }

        result.toString()
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        final result = new StringBuilder()

        for( String it : outputFiles ) {
            result
                    .append(copyMany(it, targetDir))
                    .append('\n')
        }

        return result.toString()
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

    String copyMany(String local, Path target) {
        "IFS=\$'\\n'; for name in \$(eval \"ls -1d ${Escape.path(local)}\" 2>/dev/null);do gsutil -m -q cp -R \$name ${Escape.uriPath(target)}/\$name; done; unset IFS"
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
