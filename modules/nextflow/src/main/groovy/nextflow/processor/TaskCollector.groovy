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
/**
 * Collect tasks and submit them as aggregate jobs to the underlying
 * executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class TaskCollector {

    /**
     * The set of directives which are used by the aggregate job.
     */
    private static final List<String> SUBMIT_DIRECTIVES = [
            'accelerator',
            'arch',
            'clusterOptions',
            'cpus',
            'disk',
            'machineType',
            'memory',
            'queue',
            'resourceLabels',
            'resourceLimits',
            'time',
            // only needed for container-native executors and/or Fusion
            'container',
            'containerOptions',
    ]

    private Executor executor

    private int size

    private Lock sync = new ReentrantLock()

    private List<TaskRun> batch

    private boolean closed = false

    TaskCollector(Executor executor, int size) {
        this.executor = executor
        this.size = size
        this.batch = new ArrayList<>(size)
    }

    /**
     * Add a task to the current batch, and submit the batch when it
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

            // add task to the batch
            batch.add(task)

            // submit aggregate job when it is ready
            if( batch.size() == size ) {
                executor.submit(createAggregateTask(batch))
                batch = new ArrayList<>(size)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial aggregate job.
     */
    void close() {
        sync.lock()
        try {
            if( batch.size() == 1 ) {
                executor.submit(batch.first())
            }
            else if( batch.size() > 0 ) {
                executor.submit(createAggregateTask(batch))
                batch = null
            }
            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    protected Executor getExecutor() {
        return executor
    }

    /**
     * Create the task config for an aggregate job
     * which only includes directives related to submission.
     *
     * @param processor
     */
    protected TaskConfig createAggregateConfig(TaskProcessor processor) {
        final configProps = new HashMap<String,Object>(SUBMIT_DIRECTIVES.size())
        for( final key : SUBMIT_DIRECTIVES ) {
            final value = processor.config.get(key)
            if( value != null )
                configProps[key] = value
        }
        return new TaskConfig(configProps)
    }

    /**
     * Create the task run for an aggregate job.
     *
     * @param tasks
     */
    abstract protected TaskRun createAggregateTask(List<TaskRun> tasks)

}
