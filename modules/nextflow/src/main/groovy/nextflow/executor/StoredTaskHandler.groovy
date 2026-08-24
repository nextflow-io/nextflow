/*
 * Copyright 2013-2026, Seqera Labs
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

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.TraceRecord

/**
 * Implements a {@link TaskHandler} instance for nextflow stored task ie.
 * tasks whose execution is skipped due the use of the `storeDir` directive.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class StoredTaskHandler extends TaskHandler {

    StoredTaskHandler(TaskRun task) {
        super(task)
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
    protected void killTask() {
        throw new UnsupportedOperationException()
    }

    @Override
    void submit() {
        throw new UnsupportedOperationException()
    }

    @Override
    String getStatusString() {
        "STORED"
    }

    /**
     * @return Stored tasks are not supposed to have a trace record, therefore returns {@code null}
     */
    @Override
    TraceRecord getTraceRecord() {
        return null
    }
}
