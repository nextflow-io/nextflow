/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.processor.array

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.executor.TaskArrayAware
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.processor.TaskStatus
import nextflow.util.CacheHelper
import nextflow.util.Escape

/**
 * Hold an collection of {@link TaskHandler} executed via a sole Job array submission
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Deprecated
class ArrayTaskHandler extends TaskHandler {

    private Path launcherPath

    private String launcherId

    ArrayTaskHandler(TaskRun task) {
        super(task)
        assert task.arrayTasks, "Invalid task array record"
    }

    @Override
    boolean checkIfRunning() {
        return false
    }

    @Override
    boolean checkIfCompleted() {
        return false
    }

    @Override
    void kill() {

    }

    protected List<TaskRun> getArray() { getTask().getArrayTasks() }

    @Override
    Path prepareLauncher() {

        final allLaunchers = new ArrayList<Path>(array.size())
        final script = new StringBuilder()
        // create the *array* listing all task launcher scripts
        script.append('#!/bin/bash\n')
        script.append('tasks=()')
        for( TaskRun handler : array ) {
            final l0 = handler.prepareLauncher()
            allLaunchers.add(l0)
            script.append("tasks+=(${Escape.path(l0)})")
        }
        // the launcher entry point
        final index = array.first().taskArrayIndexVariable()
        if( index ) {
            script << '# execute the i-th task\n'
            script << '/bin/bash ${tasks[' << index << ']}\n'
        }
        else {
            script << '# execute all tasks\n'
            script << 'for t in ${tasks[@]}; do /bin/bash $t; done\n'
        }
        // compute unique id & work dir for the array launcher
        final baseDir = task.processor.executor.workDir
        final hash = CacheHelper.hasher(allLaunchers).hash()
        final result = FileHelper.getWorkFolder(baseDir, hash).resolve(TaskRun.CMD_RUN)
        // save to  file
        result.text = script.toString()
        return launcherPath = result
    }

    @Override
    void submit() {
        if( !isFull() ) {
            // do nothing
            return
        }
        // do real submission
        task.submitTaskArray(launcherPath, this)

        // mark as submitted
        for( TaskArrayAware handler : array )
            handler.status = TaskStatus.SUBMITTED
    }

}
