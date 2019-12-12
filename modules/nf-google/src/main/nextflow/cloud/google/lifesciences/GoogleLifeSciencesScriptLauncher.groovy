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

import groovy.transform.CompileStatic
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun

/**
 * Implements BASH launcher script for Google Pipelines.
 *
 * @author Paolo Di Tommaso
 * @author Ólafur Haukur Flygenring <olafurh@wuxinextcode.com>
 */
@CompileStatic
class GoogleLifeSciencesScriptLauncher extends BashWrapperBuilder {

    GoogleLifeSciencesScriptLauncher(TaskBean bean, GoogleLifeSciencesTaskHandler handler) {
        super(bean, new GoogleLifeSciencesFileCopyStrategy(bean, handler))
        // enable the copying of output file to the GS work dir
        scratch = bean.workDir.toString()
        // include task script as an input to force its staging in the container work directory
        bean.inputFiles[TaskRun.CMD_SCRIPT] = bean.workDir.resolve(TaskRun.CMD_SCRIPT)
        // include task stdin file
        if (bean.input != null) {
            bean.inputFiles[TaskRun.CMD_INFILE] = bean.workDir.resolve(TaskRun.CMD_INFILE)
        }
    }

    @Override
    protected String getScratchDirectoryCommand() {
        "NXF_SCRATCH=$scratch"
    }

    @Override
    protected String getUnstageControls() {
        def result = super.getUnstageControls()
        result += copyFile(TaskRun.CMD_EXIT, workDir.resolve(TaskRun.CMD_EXIT)) + ' || true\n'
        return result
    }

    @Override
    String getCleanupCmd(String scratch) { null }

    @Override
    protected String getStageCommand() { null }

    @Override
    protected String getUnstageCommand() { null }

    @Override
    String touchFile(Path file) {
        return null
    }
}