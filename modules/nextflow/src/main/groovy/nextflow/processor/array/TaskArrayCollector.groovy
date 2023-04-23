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

import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import nextflow.processor.TaskRun

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskArrayCollector {

    final int length

    private List<TaskRun> batch

    private Lock sync = new ReentrantLock()

    TaskArrayCollector(int len) {
        this.length = len
        this.batch = new ArrayList<>(len)
    }

    /**
     * Add a task to an array batch
     *
     * @param task The task to be added
     * @return An list of {@link TaskRun} object when the batch is complete or {@code null} otherwise
     */
    List<TaskRun> batch(TaskRun task) {
        sync.lock()
        try {
            if( batch==null )
                throw new IllegalStateException("Task array collector has been already flushed")
            batch.add(task)
            if( batch.size()==length ) {
                final result = batch
                this.batch = new ArrayList<>(length)
                return result
            }
            return null
        }
        finally {
            sync.unlock()
        }
    }

    List<TaskRun> flush() {
        sync.lock()
        try {
            final result = this.batch
            this.batch = null
            return result
        }
        finally {
            sync.unlock()
        }
    }

}
