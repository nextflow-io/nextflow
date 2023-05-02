/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.processor

import java.nio.file.Files

import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.file.FileHelper
import nextflow.util.CacheHelper

/**
 * Models a task group, which executes a set of tasks sequentially
 * on the same node.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class TaskGroup extends TaskRun {

    List<TaskRun> children

    TaskGroup(List<TaskRun> tasks, Executor executor, String script) {
        this.children = tasks
        this.script = script

        // use first task by default
        final first = tasks.first()

        this.id = first.id
        this.index = first.index
        this.processor = first.processor
        this.type = first.type
        this.config = processor.config.createTaskConfig()
        this.context = new TaskContext(processor)

        // compute hash and work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)

        Files.createDirectories(workDir)

        this.hash = hash
        this.workDir = workDir
    }

    void finalize() {
        for( TaskRun task : children ) {
            task.exitStatus = exitStatus
            task.error = error
            task.stdout = task.workDir.resolve(TaskRun.CMD_OUTFILE)
            task.stderr = task.workDir.resolve(TaskRun.CMD_ERRFILE)
        }
    }

}
