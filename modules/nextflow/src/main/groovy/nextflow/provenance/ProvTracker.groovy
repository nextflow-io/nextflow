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

package nextflow.provenance

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ProvTracker {

    static class NullMessage { }

    private Map<Integer,TaskRun> messages = new ConcurrentHashMap<>()

    List<Object> beforeRun(TaskRun task, List messages) {
        // find the upstream tasks id
        findUpstreamTasks(task, messages)
        // the second entry of messages list represent the task inputs list
        // apply the de-normalization before returning it
        return denormalizeMessages(messages)
    }

    protected void findUpstreamTasks(TaskRun task, List messages) {
        // find upstream tasks and restore nulls
        final result = new HashSet<TaskId>()
        for( Object msg : messages ) {
            if( msg==null )
                throw new IllegalArgumentException("Message cannot be a null object")
            final msgId = System.identityHashCode(msg)
            result.addAll(findUpstreamTasks0(msgId,result))
        }
        // finally bind the result to the task record
        task.upstreamTasks = result
    }

    protected Set<TaskId> findUpstreamTasks0(final int msgId, Set<TaskId> upstream) {
        final task = messages.get(msgId)
        if( task==null ) {
            return upstream
        }
        if( task ) {
            upstream.add(task.id)
            return upstream
        }
        return upstream
    }

    protected Object denormalizeMessage(Object msg) {
        return msg !instanceof NullMessage ? msg : null
    }

    protected List<Object> denormalizeMessages(List messages) {
        return messages.collect(it-> denormalizeMessage(it))
    }

    protected Object normalizeMessage(Object message) {
        // map a "null" value into an  instance of "NullMessage"
        // because it's needed the object identity to track the message flow
        return message!=null ? message : new NullMessage()
    }

    void bindOutput(TaskRun task, DataflowWriteChannel ch, Object msg) {
        final value = normalizeMessage(msg)
        // map the message with the run where it has been output
        messages.put(System.identityHashCode(value), task)
        // now emit the value
        ch.bind(value)
    }

}
