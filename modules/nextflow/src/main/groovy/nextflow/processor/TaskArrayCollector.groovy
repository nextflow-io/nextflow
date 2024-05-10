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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.executor.TaskArrayExecutor
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.Escape
/**
 * Collect tasks and submit them as job arrays to the underlying
 * executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskArrayCollector extends TaskCollector {

    private TaskProcessor processor

    TaskArrayCollector(TaskProcessor processor, Executor executor, int arraySize) {
        super(executor, arraySize)
        if( executor !instanceof TaskArrayExecutor )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support job arrays")
        this.processor = processor
    }

    private TaskArrayExecutor arrayExecutor() {
        (TaskArrayExecutor)executor
    }

    /**
     * Create the task run for a job array.
     *
     * @param tasks
     */
    @Override
    protected TaskRun createAggregateTask(List<TaskRun> tasks) {
        // prepare child job launcher scripts
        final handlers = tasks.collect( t -> executor.createTaskHandler(t) )
        for( TaskHandler handler : handlers ) {
            handler.prepareLauncher()
        }

        // create work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)
        Files.createDirectories(workDir)

        // create wrapper script
        final script = createArrayTaskScript(handlers)
        log.debug "Creating task array run >> $workDir\n$script"

        // create job array
        final first = tasks.min( t -> t.index )
        final taskArray = new TaskArrayRun(
            id: first.id,
            index: first.index,
            processor: processor,
            type: processor.taskBody.type,
            config: createAggregateConfig(processor),
            context: new TaskContext(processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: handlers
        )
        taskArray.config.context = taskArray.context
        taskArray.config.process = taskArray.processor.name
        taskArray.config.executor = taskArray.processor.executor.name

        return taskArray
    }

    /**
     * Create the wrapper script for a job array.
     *
     * @param array
     */
    protected String createArrayTaskScript(List<TaskHandler> array) {
        final workDirs = array.collect( h -> Escape.path(executor.getChildWorkDir(h)) )
        """
        array=( ${workDirs.join(' ')} )
        export nxf_array_task_dir=${getArrayIndexRef()}
        ${executor.getChildLaunchCommand('$nxf_array_task_dir')}
        """.stripIndent().leftTrim()
    }

    protected String getArrayIndexRef() {
        final name = arrayExecutor().getArrayIndexName()
        final start = arrayExecutor().getArrayIndexStart()
        final index = start > 0 ? "${name} - ${start}" : name
        return '${array[' + index + ']}'
    }

}
