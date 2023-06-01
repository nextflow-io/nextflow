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
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.executor.TaskArrayAware
import nextflow.file.FileHelper
import nextflow.util.CacheHelper
import nextflow.util.Escape

/**
 * Task monitor that batches tasks and submits them as array jobs
 * to an underlying task monitor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskArrayCollector {

    private TaskArrayAware executor

    private int arraySize

    private Lock sync = new ReentrantLock()

    private List<TaskHandler> array

    private boolean closed = false

    TaskArrayCollector(Executor executor, int arraySize) {
        if( executor !instanceof TaskArrayAware )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support array jobs")

        this.executor = (TaskArrayAware)executor
        this.arraySize = arraySize
        this.array = new ArrayList<>(arraySize)
    }

    /**
     * Add a task to the current array, and submit the array when it
     * reaches the desired size.
     *
     * @param task
     */
    void collect(TaskRun task) {
        sync.lock()

        try {
            // submit task directly if the collector is closed
            // or if the task is retried (since it might have dynamic resources)
            if( closed || task.config.getAttempt() > 1 ) {
                executor.submit(task)
                return
            }

            // create task handler
            final handler = executor.createTaskHandler(task)

            // add task to the array
            array << handler

            // submit array job when it is ready
            if( array.size() == arraySize ) {
                submit0(array)
                array = new ArrayList<>(arraySize)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial array job.
     */
    void close() {
        sync.lock()

        try {
            if( array.size() > 0 )
                submit0(array)

            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    protected void submit0(List<TaskHandler> array) {
        // prepare child job launcher scripts
        for( TaskHandler handler : array )
            handler.prepareLauncher()

        // submit array job
        executor.submit(createTaskArray(array))
    }

    /**
     * Create the task run for an array job.
     *
     * @param array
     */
    protected TaskRun createTaskArray(List<TaskHandler> array) {
        final tasks = array.collect( h -> h.task )
        final first = tasks.first()

        // create work directory
        final hash = CacheHelper.hasher( tasks.collect( t -> t.getHash().asLong() ) ).hash()
        final workDir = FileHelper.getWorkFolder(executor.getWorkDir(), hash)

        Files.createDirectories(workDir)

        // create wrapper script
        final script = createTaskArrayScript(array)

        // create task handler
        return new TaskArray(
            id: first.id,
            index: first.index,
            processor: first.processor,
            type: first.type,
            config: first.processor.config.createTaskConfig(),
            context: new TaskContext(first.processor),
            hash: hash,
            workDir: workDir,
            script: script,
            children: array
        )
    }

    /**
     * Create the wrapper script for an array job.
     *
     * @param array
     */
    protected String createTaskArrayScript(List<TaskHandler> array) {
        // get work directory and launch command for each task
        final workDirs = array.collect( h -> h.getWorkDir() )
        final args = array.first().getLaunchCommand().toArray() as String[]
        final cmd = Escape.cli(args).replaceAll(workDirs.first(), '\\${task_dir}')

        // create wrapper script
        final arrayIndexName = executor.getArrayIndexName()
        final builder = new StringBuilder()
            << "array=( ${workDirs.collect( p -> Escape.path(p) ).join(' ')} )\n"
            << "export task_dir=\${array[${arrayIndexName}]}\n"
            << cmd << '\n'
        return builder.toString()
    }

}
