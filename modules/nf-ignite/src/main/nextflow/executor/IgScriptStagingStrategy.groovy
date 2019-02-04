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

package nextflow.executor

import java.nio.file.Files
import java.nio.file.Path

import nextflow.processor.TaskRun

/**
 * Implements a {@link StagingStrategy} interface extending the behavior of {@link IgFileStagingStrategy}
 * so that it copies also the task control meta-files.
 *
 * Moreover it extends the {@link SimpleFileCopyStrategy} class in such way that {@link ScriptFileCopyStrategy#getStageInputFilesScript}
 * and {@linl ScriptFileCopyStrategy#getUnstageOutputFilesScript} returns an empty string (because it is not required to
 * manage stage/unstage at wrapper script level)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgScriptStagingStrategy extends IgFileStagingStrategy implements ScriptFileCopyStrategy {

    @Delegate(interfaces=false)
    SimpleFileCopyStrategy delegate

    /**
     * {@inheritDoc}
     */
    @Override
    String getBeforeStartScript() {
        // do not use delegate `SimpleFileCopyStrategy#getBeforeStartScript` because
        // it uses the Global.session which is not accessible from Ignite remote nodes
        return null
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void stage() {
        delegate = new SimpleFileCopyStrategy(task)
        super.stage()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    void unstage() {
        super.unstage()
        // copy the 'exit' file and 'output' file
        copyFromScratchToWorkDir(TaskRun.CMD_EXIT)
        copyFromScratchToWorkDir(TaskRun.CMD_OUTFILE)
        copyFromScratchToWorkDir(TaskRun.CMD_ERRFILE, true)
        copyFromScratchToWorkDir(TaskRun.CMD_TRACE, true)
        copyFromScratchToWorkDir(TaskRun.CMD_LOG, true)
    }


    private void copyFromScratchToWorkDir( String name, boolean ignoreError = false ) {
        try {
            Files.copy(localWorkDir.resolve(name), task.workDir.resolve(name))
        }
        catch( Exception e ) {
            if( !ignoreError )
                log?.debug "Unable to copy file: '$name' from: '$localWorkDir' to: '$task.workDir'"
        }
    }

    /**
     * Turn off file staging from the wrapper script returning a null string
     *
     * @return A {@code null} string
     */
    @Override
    String getStageInputFilesScript(Map<String,Path> inputFiles) {
        return null
    }

    /**
     * Turn off remote file download
     *
     * @param inputFiles
     * @return
     */
    @Override
    Map<String,Path> resolveForeignFiles(Map<String,Path> inputFiles) {
        return inputFiles
    }

    /**
     * Turn off file unstage from the wrapper script returning a null string
     *
     * @return A {@code null} string
     */
    @Override
    String getUnstageOutputFilesScript(List<String> outputFiles, Path targetDir) {
        return null
    }

}
