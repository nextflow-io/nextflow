/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.processor

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TaskTracker extends DataflowEventAdapter {

    static class NullMessage { }

    @Canonical
    static class RunEntry {
        int runId
        /** Both process and operator runs have associated a {@link DataflowProcessor} */
        DataflowProcessor processor
        /** The corresponding {@link TaskStartParams} in the case process run */
        TaskStartParams task
        /** The identities of the upstream tasks */
        Set<TaskId> upstreamIds
    }

    private static final Map<Integer,RunEntry> messages = new ConcurrentHashMap<>()

    private static final ThreadLocal<RunEntry> entry = new ThreadLocal<>()

    private static final AtomicInteger generator = new AtomicInteger()

    private static final Map<TaskId,RunEntry> runs = new ConcurrentHashMap<>()

    static int len() {
        return messages.size()
    }

    static void clean() {
        messages.clear()
        generator.set(0)
    }

    @Override
    List<Object> beforeRun(DataflowProcessor processor, List<Object> messages) {
        assert messages.size() == 2, "Unexpected processes input messages size=${messages.size()}"
        assert messages[0] instanceof TaskStartParams
        assert messages[1] instanceof List

        final task = getTaskParams(messages)
        final entry = runs.get(task.id)
        if( entry==null ) {
            // nothing to do
            // this is expected for the first task that has not yet emitted any output message
            return messages
        }
        // find the upstream tasks id
        final upstream = findUpstreamTasks(messages)
        // the second entry of messages list represent the task inputs list
        // apply the de-normalization before returning it
        final inputs = messages[1] as List
        messages[1] = denormalizeMessages(inputs)
        log.debug "Input messages ${inputs.collect(it-> str(it)).join(', ')}\n^^ Upstream tasks: ${upstream.join(',')}"
        return messages
    }

    protected String str(Object o) {
        return o!=null ? "${o.class.name}@${System.identityHashCode(o)}" : null
    }

    protected Set<TaskId> findUpstreamTasks(final int msgId, Set<TaskId> upstream = new HashSet<>(10)) {
        final entry = messages.get(msgId)
        if( entry==null ) {
            return upstream
        }
        if( entry.task ) {
            upstream.add(entry.task.id)
            return upstream
        }
        return upstream
    }

    protected Set<TaskId> findUpstreamTasks(List messages) {
        // find upstream tasks and restore nulls
        final result = new HashSet<TaskId>()
        for( Object msg : messages ) {
            if( msg==null )
                throw new IllegalArgumentException("Message cannot be a null object")
            final msgId = System.identityHashCode(msg)
            result.addAll(findUpstreamTasks(msgId))
        }
        return result
    }

    @Override
    Object messageSentOut(DataflowProcessor processor, DataflowWriteChannel<Object> channel, int index, Object message) {
        // normalize the message
        final msg = normalizeMessage(message)
        // determine the current run entry
        final entry = getOrCreateEntry(processor)
        // map the message with the run entry where it has been output
        messages.put(System.identityHashCode(msg), entry)
        return msg
    }

    protected RunEntry getOrCreateEntry(DataflowProcessor processor) {
        def entry = entry.get()
        if( entry==null ) {
            entry = new RunEntry(generator.incrementAndGet(), processor)
            TaskTracker.entry.set(entry)
        }
        return entry
    }

    protected Object normalizeMessage(Object message) {
        // map a "null" value into an  instance of "NullMessage"
        // because it's needed the object identity to track the message flow
        return message!=null ? message : new NullMessage()
    }

    protected Object denormalizeMessage(Object msg) {
        return msg !instanceof NullMessage ? msg : null
    }

    protected List<Object> denormalizeMessages(List messages) {
        return messages.collect(it-> denormalizeMessage(it))
    }

    protected TaskStartParams getTaskParams(List messages) {
        if( messages && messages[0] instanceof TaskStartParams )
            return messages[0] as TaskStartParams
        else
            throw new IllegalArgumentException("Unable to find required TaskStartParams object - offending messages=$messages")
    }

    /**
     * track all outputs send out
     */
    @Override
    void afterRun(DataflowProcessor processor, List<Object> inputs) {
        final entry = entry.get()
        assert entry!=null, "Missing run entry for process=$processor"
        try {
            // update the entry with the task param
            entry.task = getTaskParams(inputs)
            // map the run entry with the task Id
            runs.put(entry.task.id, entry)
        }
        finally {
            // cleanup the current run entry
            TaskTracker.entry.remove()
        }
    }
}
