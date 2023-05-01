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

import java.nio.file.Path
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.ArrayTaskAware
import nextflow.util.Escape
/**
 * Task monitor that batches tasks and submits them as array jobs
 * to an underlying task monitor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskCollector {

    private TaskProcessor processor

    private ArrayTaskAware executor

    private int arraySize

    private Lock sync = new ReentrantLock()

    private List<TaskRun> array

    private volatile boolean closed

    ArrayTaskCollector(TaskProcessor processor, int arraySize) {
        if( processor.executor !instanceof ArrayTaskAware )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support array jobs")
        this.processor = processor
        this.executor = (ArrayTaskAware) processor.executor
        this.arraySize = arraySize
        this.array = new ArrayList<>(arraySize)
    }

    /**
     * Add a task to the current array, and submit the array when it
     * reaches the desired size.
     *
     * @param handler
     */
    void submit(TaskRun task) {
        // submit task directly if the collector is closed
        if( closed ) {
            throw new IllegalStateException("Process is already closed")
        }

        sync.lock()
        try {
            // update the array index attribute
            task.arrayIndex = array.size()
            // add task to the array array
            array.add(task)

            // submit array job when the array is ready
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
     *
     * @param process
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

    protected void submit0(List<TaskRun> array) {
        // final submit the task array
        final task = processor.createTaskArray(array)
        task.script = createTaskArrayScript(array)
        executor.submit(task)
    }

    protected String createTaskArrayScript(List<TaskRun> array) {
        final allLaunchers = new ArrayList<Path>()
        final script = new StringBuilder()
        script.append('tasks=()').append('\n')
        for( TaskRun t0 : array ) {
            final l0 = t0.workDir.resolve(TaskRun.CMD_RUN)
            allLaunchers.add(l0)
            script.append("tasks+=(${Escape.path(l0)})").append('\n')
        }
        // the launcher entry point
        final index = executor.getArrayIndexName()
        if( index ) {
            script << '# execute the i-th task\n'
            script << '/bin/bash ${tasks[' << index << ']}\n'
        }
        else {
            script << '# execute all tasks\n'
            script << 'for t in ${tasks[@]}; do /bin/bash $t; done\n'
        }
    }

}
