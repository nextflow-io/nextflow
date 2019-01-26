/*
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

package nextflow.cloud.google.pipelines

import java.nio.file.Path
import java.nio.file.Paths

import com.google.cloud.storage.contrib.nio.CloudStoragePath
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.util.Escape

/**
 * Defines the file/script copy strategies for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@Slf4j
@CompileStatic
class GooglePipelinesFileCopyStrategy extends SimpleFileCopyStrategy {

    GooglePipelinesTaskHandler handler
    TaskBean task

    GooglePipelinesFileCopyStrategy(TaskBean bean, GooglePipelinesTaskHandler handler) {
        super(bean)
        this.handler = handler
        this.task = bean
    }

    @Override
    String getStageInputFilesScript(Map<String, Path> inputFiles) {

        def createDirectories  = []
        def stagingCommands = []

        inputFiles.each { stageName, storePath ->
            def absStorePath = storePath.toAbsolutePath()
            def storePathIsDir = absStorePath.isDirectory()
            def stagePath = Paths.get(stageName)
            def parent = stagePath.parent
            //we can't simply use the handy ToUriString function since it also URL encodes spaces and gsutil doesn't like that
            def escapedStoreUri = "${absStorePath.getFileSystem()}${Escape.path(absStorePath)}"
            def escapedStageName = Escape.path(stageName)

            //check if we need to create parent dirs for staging file since gsutil doesn't create them for us
            if(parent) {
                createDirectories << "mkdir -p $workDir/$parent".toString()
            }

            if(storePathIsDir) {
                stagingCommands << "gsutil -m -q cp -R $escapedStoreUri/ $workDir".toString()
                //check if we need to move the directory (gsutil doesn't support renaming directories on copy)
                if(parent || !absStorePath.toString().endsWith(stageName)) {
                    def newLocation = "$workDir/$escapedStageName"
                    stagingCommands << "mv $workDir/${Escape.path(absStorePath.name)} $newLocation".toString()
                }
            } else {
                stagingCommands << "gsutil -m -q cp $escapedStoreUri $workDir/$escapedStageName".toString()
            }
        }

        handler.stagingCommands.addAll(createDirectories)
        handler.stagingCommands.addAll(stagingCommands)

        // copy the remoteBinDir if it ia defined
        if(handler.pipelineConfiguration.remoteBinDir) {
            def createRemoteBinDir = "mkdir -p $workDir/nextflow-bin".toString()
            def remoteBinCopy = "gsutil -m -q cp -P -r ${handler.pipelineConfiguration.remoteBinDir.toUriString()}/* $workDir/nextflow-bin".toString()
            handler.stagingCommands << createRemoteBinDir << remoteBinCopy
        }

        // Insert this comment into the task run script to note that the staging is done differently
        "# Google pipeline staging is done in a pipeline action step that is run prior to the main pipeline action"
    }

    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {

        def unstagingCommands = outputFiles.collect {
            def escaped = Escape.path(it)
            "gsutil -m -q cp -R $escaped ${workDir.toUriString()} || true".toString()
        }

        log.debug "[GPAPI] Constructed the following file copy staging commands: $unstagingCommands"

        handler.unstagingCommands.addAll(unstagingCommands)

        //Insert this comment into the task run script to note that the unstaging is done differently
        ": # Google pipeline unstaging is done in a pipeline action step that is run after the main pipeline action"
    }

    //Although it seems like this is used as a general "touch" mechanism it's actually just used to create a .begin in the working directory
    @Override
    String touchFile(Path file) {
        if (file in CloudStoragePath) {
            handler.stagingCommands << "echo start | gsutil -q cp  -c - ${file.toUriString()} || true".toString()
            "# Google pipeline touchFile is done in a pipeline action step that is run prior to the main pipeline action"
        } else
            return super.touchFile(file)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    String getEnvScript(Map environment, boolean container) {
        if( container )
            throw new IllegalArgumentException("Parameter `wrapHandler` not supported by ${this.class.simpleName}")

        final result = new StringBuilder()
        final copy = environment ? new HashMap<String,String>(environment) : Collections.<String,String>emptyMap()
        // remove any external PATH
        if( copy.containsKey('PATH') )
            copy.remove('PATH')
        // when a remote bin directory is provide managed it properly
        if( handler.pipelineConfiguration.remoteBinDir ) {
            result << "chmod +x $workDir/nextflow-bin/*\n"
            result << "export PATH=$workDir/nextflow-bin:\$PATH\n"
        }
        // finally render the environment
        final envSnippet = super.getEnvScript(copy,false)
        if( envSnippet )
            result << envSnippet
        return result.toString()
    }
}