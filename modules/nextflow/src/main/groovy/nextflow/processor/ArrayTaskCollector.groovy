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

import java.util.concurrent.ConcurrentLinkedQueue

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.ArrayTaskAware
import nextflow.executor.Executor

/**
 * Task monitor that batches tasks and submits them as array jobs
 * to an underlying task monitor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskCollector {

    private Executor executor

    private TaskMonitor monitor

    private int arraySize

    private volatile Queue<TaskHandler> queue = new ConcurrentLinkedQueue<>()

    private volatile boolean closed = false

    ArrayTaskCollector(Executor executor, int arraySize) {
        if( executor !instanceof ArrayTaskAware )
            throw new IllegalArgumentException("Executor '${executor.name}' does not support array jobs")

        this.executor = executor
        this.monitor = executor.monitor
        this.arraySize = arraySize
    }

    /**
     * Add a task to the queue, and submit an array job when the
     * queue reaches the desired size.
     *
     * @param handler
     */
    synchronized void submit(TaskRun task) {
        // submit task directly if the collector is closed
        if( closed ) {
            executor.submit(task)
            return
        }

        // create task handler
        final handler = executor.createTaskHandler(task)

        // add task to the array queue
        queue.add(handler)

        // submit array job when a batch is ready
        if( queue.size() >= arraySize )
            submit0(queue, arraySize)
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial array job.
     *
     * @param process
     */
    void close() {
        if( queue.size() > 0 )
            submit0(queue, queue.size())

        closed = true
    }

    synchronized protected void submit0( Queue<TaskHandler> queue, int size ) {
        // remove tasks from the queue
        def array = new ArrayList<TaskHandler>(size)
        def iter = queue.iterator()

        while( iter.hasNext() && array.size() < size ) {
            array << iter.next()
            iter.remove()
        }

        if( array.size() < size )
            log.warn "Array job expected ${size} tasks but received only ${array.size()}"

        // create submitter for array job
        ((ArrayTaskAware)executor).createArrayTaskSubmitter(array)

        // submit each task to the underlying monitor
        for( TaskHandler handler : array )
            monitor.schedule(handler)
    }

}
