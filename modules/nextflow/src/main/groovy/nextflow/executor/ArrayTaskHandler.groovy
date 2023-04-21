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
 *
 */

package nextflow.executor


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.processor.TaskHandler
import nextflow.processor.TaskStatus
import nextflow.trace.TraceRecord
/**
 * Generic handler for an array task, which simply launches
 * each task and waits for it to finish.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ArrayTaskHandler extends TaskHandler {

    final List<TaskHandler> array

    ArrayTaskHandler(List<TaskHandler> array) {
        // use the first task to provide common properties (e.g. configuration)
        super(array.first().getTask())

        this.array = array
    }

    @Override
    void submit() {
        for( TaskHandler handler : array )
            handler.submit()

        status = TaskStatus.SUBMITTED
    }

    @Override
    boolean checkIfRunning() {
        for( TaskHandler handler : array )
            if( !handler.checkIfRunning() )
                return false

        status = TaskStatus.RUNNING
        return true
    }

    @Override
    boolean checkIfCompleted() {
        for( TaskHandler handler : array )
            if( !handler.checkIfCompleted() )
                return false

        status = TaskStatus.COMPLETED
        return true
    }

    @Override
    void kill() {
        for( TaskHandler handler : array )
            handler.kill()
    }

    @Override
    protected StringBuilder toStringBuilder( StringBuilder builder ) {
        builder << '\n'
        for( TaskHandler handler : array )
            builder << '    ' << handler.toString() << '\n'

        return builder
    }

    @Override
    String getStatusString() {
        if( array.any { h -> h.task.failed } ) return 'FAILED'
        if( array.any { h -> h.task.aborted } ) return 'ABORTED'
        return this.status.toString()
    }

    @Override
    TraceRecord safeTraceRecord() {
        throw new UnsupportedOperationException()
    }

    @Override
    TraceRecord getTraceRecord() {
        throw new UnsupportedOperationException()
    }

    @Override
    boolean canForkProcess() {
        final max = task.processor.maxForks
        return !max
            ? true
            : task.processor.forksCount + array.size() <= max
    }

    @Override
    void incProcessForks() {
        task.processor.forksCount?.add(array.size())
    }

    @Override
    void decProcessForks() {
        task.processor.forksCount?.add(-array.size())
    }

}
