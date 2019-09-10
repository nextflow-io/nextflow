/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord

/**
 * Implements a cached {@link TaskHandler}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CachedTaskHandler extends TaskHandler {

    private TraceRecord trace

    CachedTaskHandler(TaskRun task, TraceRecord trace) {
        super(task)
        this.trace = trace
        this.trace.setCached(true)
        // override task id with the current one
        // https://github.com/nextflow-io/nextflow/issues/1301
        this.trace.put('task_id', task.id)
    }

    @Override
    boolean checkIfRunning() {
        return false
    }

    @Override
    boolean checkIfCompleted() {
        return true
    }

    @Override
    void kill() {
        throw new UnsupportedOperationException()
    }

    @Override
    void submit() {
        throw new UnsupportedOperationException()
    }

    @Override
    String getStatusString() {
        "CACHED"
    }

    @Override
    TraceRecord getTraceRecord() {
        return trace
    }

}
