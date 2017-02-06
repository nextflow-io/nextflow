/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor

import java.nio.file.Files

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

    @Delegate
    ScriptFileCopyStrategy delegate

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
    String getStageInputFilesScript() {
        return null
    }

    /**
     * Turn off file unstage from the wrapper script returning a null string
     *
     * @return A {@code null} string
     */
    @Override
    String getUnstageOutputFilesScript() {
        return null
    }

}
