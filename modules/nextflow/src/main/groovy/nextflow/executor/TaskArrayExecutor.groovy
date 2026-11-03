/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.FileSystems
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.fusion.FusionHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
/**
 * Interface for executors that support job arrays.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
// declared as a trait to keep binary compatibility with plugins compiled with Groovy 4
trait TaskArrayExecutor {

    abstract String getName()

    abstract Path getWorkDir()

    abstract void submit( TaskRun task )

    abstract TaskHandler createTaskHandler(TaskRun task)

    abstract boolean isFusionEnabled()

    /**
     * Get the environment variable name that provides the array index of a task.
     */
    abstract String getArrayIndexName()

    /**
     * Get the start of the job array index range.
     */
    abstract int getArrayIndexStart()

    /**
     * Get the name of a child job based on the array job name
     * and child index.
     */
    abstract String getArrayTaskId(String jobId, int index)

    boolean isWorkDirDefaultFS() {
        getWorkDir().fileSystem== FileSystems.default
    }

    /**
     * Get a {@link TaskHandler} work directory for the task array resolution
     *
     * @param handler
     * @return
     */
    String getArrayWorkDir(TaskHandler handler) {
        return isFusionEnabled()
            ? FusionHelper.toContainerMount(handler.task.workDir).toString()
            : handler.task.workDir.toUriString()
    }

    String getArrayLaunchCommand(String taskDir) {
        if( isFusionEnabled() ) {
            return "bash ${taskDir}/${TaskRun.CMD_RUN}"
        }
        else if( isWorkDirDefaultFS() ) {
            return "bash ${taskDir}/${TaskRun.CMD_RUN} 2>&1 > ${taskDir}/${TaskRun.CMD_LOG}"
        }
        else {
            throw new IllegalStateException("Executor ${getName()} does not support array jobs")
        }
    }
}
