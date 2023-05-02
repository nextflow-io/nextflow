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

import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.executor.TaskArrayAware
import nextflow.executor.TaskArraySubmitter

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

    private TaskMonitor monitor

    private int arraySize

    private Lock sync = new ReentrantLock()

    private List<TaskHandler> array

    private boolean closed = false

    TaskArrayCollector(Executor executor, int arraySize) {
        if( executor !instanceof TaskArrayAware )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support array jobs")

        this.executor = (TaskArrayAware)executor
        this.monitor = executor.monitor
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
            if( closed ) {
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
        // create submitter for array job
        final arraySubmitter = new TaskArraySubmitter(array, executor)

        // submit each task to the underlying monitor
        // each task will defer to the array job during submission
        for( TaskHandler handler : array ) {
            handler.arraySubmitter = arraySubmitter
            monitor.schedule(handler)
        }
    }

}
