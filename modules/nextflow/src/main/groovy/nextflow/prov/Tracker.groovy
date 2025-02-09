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

package nextflow.prov

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.extension.op.Op
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
/**
 * Implement provenance tracking logic for tasks and operators
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class Tracker {

    static @Canonical class Msg {
        final Object value

        String toString() {
            "Msg[id=${System.identityHashCode(this)}; value=${value}]"
        }

        static Msg of(Object o) {
            new Msg(o)
        }
    }

    /**
     * Associate an output value with the corresponding task run that emitted it
     */
    private Map<Integer,ProvLink> messages = new ConcurrentHashMap<>()

    List<Object> receiveInputs(TaskRun task, List inputs) {
        // find the upstream tasks id
        findUpstreamTasks(task, inputs)
        // log for debugging purposes
        logInputs(task, inputs)
        // the second entry of messages list represent the run inputs list
        // apply the de-normalization before returning it
        return Op.unwrap(inputs)
    }

    private logInputs(TaskRun task, List inputs) {
        if( log.isTraceEnabled() ) {
            def msg = "Task input"
            msg += "\n - task id : ${task.id} "
            msg += "\n - name    : '${task.name}'"
            msg += "\n - upstream: ${task.upstreamTasks*.value}"
            for( Object it : inputs ) {
                msg += "\n=> ${dumpInput(it)}"
            }
            log.trace(msg)
        }
    }

    private logInputs(OperatorRun run, List inputs) {
        if( log.isTraceEnabled() ) {
            def msg = "Operator input"
            msg += "\n - run id: ${System.identityHashCode(run)} "
            for( Object it : inputs ) {
                msg += "\n=> ${dumpInput(it)}"
            }
            log.trace(msg)
        }
    }

    private String dumpInput(Object it) {
        if( it==null )
            return "null"
        if( it instanceof Msg )
            return it.toString()
        else
            return "${it.getClass().getName()}[id=${System.identityHashCode(it)}; value=${it.toString()}]"
    }

    List<Object> receiveInputs(OperatorRun run, List inputs) {
        // find the upstream tasks id
        run.inputIds.addAll(inputs.collect(msg-> System.identityHashCode(msg)))
        // log for debugging purposes
        logInputs(run, inputs)
        // the second entry of messages list represent the task inputs list
        // apply the de-normalization before returning it
        return Op.unwrap(inputs)
    }

    protected void findUpstreamTasks(TaskRun task, List messages) {
        // find upstream tasks and restore nulls
        final result = new HashSet<TaskId>()
        for( Object msg : messages ) {
            if( msg==null )
                throw new IllegalArgumentException("Message cannot be a null object")
            if( msg !instanceof Msg )
                continue
            final msgId = System.identityHashCode(msg)
            result.addAll(findUpstreamTasks0(msgId,result))
        }
        // finally bind the result to the task record
        task.upstreamTasks = result
    }

    protected Set<TaskId> findUpstreamTasks0(final int msgId, Set<TaskId> upstream) {
        final run = messages.get(msgId)
        if( run instanceof TaskRun ) {
            upstream.add(run.id)
            return upstream
        }
        if( run instanceof OperatorRun ) {
            for( Integer it : run.inputIds ) {
                if( it!=msgId ) {
                    findUpstreamTasks0(it, upstream)
                }
                else {
                    log.trace "Skip duplicate provenance message id=${msgId}"
                }
            }
        }
        return upstream
    }

    Msg bindOutput(ProvLink run, DataflowWriteChannel channel, Object out) {
        assert run!=null, "Argument 'run' cannot be null"
        assert channel!=null, "Argument 'channel' cannot be null"

        final msg = Op.wrap(out)
        logOutput(run, msg)
        // map the message with the run where it has been output
        messages.put(System.identityHashCode(msg), run)
        // now emit the value
        channel.bind(msg)
        return msg
    }

    private void logOutput(ProvLink run, Msg msg) {
        if( log.isTraceEnabled() ) {
            String str
            if( run instanceof OperatorRun ) {
                str = "Operator output"
                str += "\n - run id: ${System.identityHashCode(run)}"
            }
            else if( run instanceof TaskRun ) {
                str = "Task output"
                str += "\n - task id : ${run.id}"
                str += "\n - name    : '${run.name}'"
            }
            else
                throw new IllegalArgumentException("Unknown run type: ${run}")
            str += "\n<= ${msg}"
            log.trace(str)
        }
    }

}
