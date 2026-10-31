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

package nextflow.dataflow.engine

import java.util.concurrent.atomic.AtomicInteger
import java.util.function.BooleanSupplier

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import groovyx.gpars.dataflow.DataflowStreamReadAdapter

/**
 * The backlog of a process:
 * - preparing: input sets taken by the process and not yet scheduled (or discarded)
 * - pending: tasks in the executor pending queue, until accepted for submission
 * - ready: complete input sets queued at the process inputs
 *
 * The process operator pauses when {@code preparing + pending} reaches the high-water mark,
 * and resumes once it drops to the low-water mark. Sources feeding the process are paused
 * likewise on the whole backlog.
 *
 * A task counts as both preparing and pending from its scheduling ({@code schedule} is called
 * within the task preparation) until the end of its preparation ({@code afterRun}), so the
 * backlog may briefly exceed the mark by the input sets being prepared.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ProcessBacklog implements BooleanSupplier {

    final String name

    final DataflowProcessor operator

    private final Backpressure owner

    private final List<DataflowReadChannel> inputs

    private final AtomicInteger preparing = new AtomicInteger()

    private final AtomicInteger pending = new AtomicInteger()

    private volatile boolean closed

    @PackageScope
    ProcessBacklog(Backpressure owner, String name, List<DataflowReadChannel> inputs, DataflowProcessor operator) {
        this.owner = owner
        this.name = name
        this.inputs = inputs
        this.operator = operator
    }

    /**
     * The operator gate: true while the process must not take new inputs
     */
    @Override
    boolean getAsBoolean() {
        return !closed && !owner.released && preparing.get() + pending.get() >= owner.highWater
    }

    void taskStarted() {
        preparing.incrementAndGet()
    }

    void taskPrepared() {
        preparing.decrementAndGet()
        changed()
    }

    void taskPending() {
        pending.incrementAndGet()
    }

    void taskSubmitted() {
        pending.decrementAndGet()
        changed()
    }

    /**
     * Close the backlog when the process terminates, and deregister it from the backpressure
     */
    void close() {
        closed = true
        changed()
        owner.removeProcess(this)
    }

    private void changed() {
        if( preparing.get() + pending.get() <= owner.lowWater )
            operator?.resume()
        owner.wake()
    }

    boolean isPaused() { operator?.isPaused() }

    /**
     * The number of times the process operator was paused
     */
    @PackageScope
    long getPauseCount() { operator ? operator.getPauseCount() : 0 }

    int getPreparing() { preparing.get() }

    int getPending() { pending.get() }

    /**
     * Whether the whole backlog is at or above the given mark
     */
    boolean isOver(int mark) {
        if( closed )
            return false
        final busy = preparing.get() + pending.get()
        return busy >= mark || busy + ready(mark - busy) >= mark
    }

    /**
     * The number of complete input sets queued at the process, counting at most {@code max}.
     * Items waiting for an unbound value input are not counted.
     */
    int ready(int max) {
        int result = max
        boolean anyQueue = false
        for( DataflowReadChannel ch : inputs ) {
            if( ch instanceof DataflowVariable ) {
                if( !ch.isBound() )
                    return 0
                continue
            }
            anyQueue = true
            result = Math.min(result, queued(ch, result))
            if( result == 0 )
                return 0
        }
        return anyQueue ? result : 0
    }

    static private int queued(DataflowReadChannel ch, int max) {
        if( ch instanceof DataflowStreamReadAdapter )
            return ch.length(max)
        if( ch instanceof DataflowQueue )
            return Math.min(ch.length(), max)
        return 0
    }

    @Override
    String toString() {
        final hw = owner.highWater
        return "${name}: preparing=${preparing.get()}; pending=${pending.get()}; ready=${ready(hw)}${ready(hw)>=hw ? '+' : ''}; paused=${isPaused()}; closed=${closed}"
    }
}
