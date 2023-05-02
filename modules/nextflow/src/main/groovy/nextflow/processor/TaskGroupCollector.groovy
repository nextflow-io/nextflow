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
 * Collect tasks and submit them as task groups to the underlying
 * executor.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskGroupCollector {

    private Executor executor

    private int groupSize

    private Lock sync = new ReentrantLock()

    private List<TaskRun> group

    private boolean closed = false

    TaskGroupCollector(Executor executor, int groupSize) {
        this.executor = executor
        this.groupSize = groupSize
        this.group = new ArrayList<>(groupSize)
    }

    /**
     * Add a task to the current group, and submit the group when it
     * reaches the desired size.
     *
     * @param task
     */
    void collect(TaskRun task) {
        sync.lock()

        try {
            // submit task directly if the collector is closed
            if( closed )
                executor.submit(task)

            // add task to the group
            group << task

            // submit task group when it is ready
            if( group.size() == groupSize ) {
                submit0(group)
                group = new ArrayList<>(groupSize)
            }
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Close the collector, submitting any remaining tasks as a partial task group.
     */
    void close() {
        sync.lock()

        try {
            if( group.size() > 0 )
                submit0(group)

            closed = true
        }
        finally {
            sync.unlock()
        }
    }

    protected void submit0(List<TaskRun> tasks) {
        // prepare work directory for each child task
        for( TaskRun task : tasks )
            executor.createTaskHandler(task).prepareLauncher()

        // create task group
        final taskGroup = new TaskGroup(tasks, executor)

        // submit task group to the underlying executor
        executor.submit(taskGroup)
    }

}
