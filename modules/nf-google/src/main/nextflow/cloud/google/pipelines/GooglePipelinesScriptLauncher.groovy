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

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

import java.nio.file.Path

/**
 * Implements BASH launcher script for Google Pipelines.
 *
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@CompileStatic
class GooglePipelinesScriptLauncher extends BashWrapperBuilder {

    private GooglePipelinesConfiguration pipelineConfiguration

    GooglePipelinesScriptLauncher(TaskBean bean, GooglePipelinesTaskHandler handler) {
        super(bean, new GooglePipelinesFileCopyStrategy(bean, handler))
        this.pipelineConfiguration = handler.pipelineConfiguration
        // enable the copying of output file to the GS work dir
        scratch = "$bean.workDir/scratch".toString()
        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // include the wrapper script as in input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_RUN] = bean.workDir.resolve(TaskRun.CMD_RUN)
        // include task stdin file
        if (bean.input != null) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }
    }

    @Override
    protected String getScratchDirectoryCommand() {
        "NXF_SCRATCH=\$(mkdir $scratch)"
    }

    @Override
    String copyFile(String name, Path target) {
        "echo 'Google Pipelines file staging/unstaging happens in pre/post actions'"
    }

}