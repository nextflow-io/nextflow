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

package nextflow.extension

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.NF
import nextflow.dag.NodeMarker
import nextflow.exception.ScriptRuntimeException
import nextflow.script.ChainableDef
import nextflow.script.ChannelOut
import nextflow.script.ComponentDef
import nextflow.script.CompositeDef
import nextflow.script.ExecutionStack
import org.codehaus.groovy.runtime.InvokerHelper
import static nextflow.util.LoggerHelper.fmtType

/**
 * Implements dataflow channel extension methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelEx {

//    static void set(ChannelOut source, Closure holder) {
//        final names = CaptureProperties.capture(holder)
//        if( names.size() > source.size() )
//            throw new IllegalArgumentException("Operation `set` expects ${names.size()} channels but only ${source.size()} are provided")
//
//        for( int i=0; i<source.size(); i++ ) {
//            final ch = source[i]
//            final nm = names[i]
//            NF.binding.setVariable(nm, ch)
//        }
//    }

    static DataflowWriteChannel dump(final DataflowWriteChannel source, Closure closure = null) {
        dump(source, Collections.emptyMap(), closure)
    }

    static DataflowWriteChannel dump(final DataflowWriteChannel source, Map opts, Closure closure = null) {
        def op = new DumpOp(opts, closure)
        if( op.isEnabled() ) {
            op.setSource(source)
            def target = op.apply()
            NodeMarker.addOperatorNode('dump', source, target)
            return target
        }
        else {
            return source
        }
    }

    /**
     * Apply a channel as the first argument to a closure and return the result.
     *
     * @param source
     * @param closure
     */
    static Object then(final DataflowWriteChannel source, Closure closure) {
        def out = closure.call(source)
        if( out instanceof DataflowWriteChannel || out instanceof ChannelOut )
            return out
        throw new ScriptRuntimeException("Pipeline closure did not return a channel or multi-channel")
    }

    static Object then(final ChannelOut source, Closure closure) {
        def out = closure.call(source)
        if( out instanceof DataflowWriteChannel || out instanceof ChannelOut )
            return out
        throw new ScriptRuntimeException("Pipeline closure did not return a channel or multi-channel")
    }

    /**
     * Creates a channel emitting the entries in the collection to which is applied
     *
     * @param values A list of values to be emitted by the resulting channel
     * @return
     */
    static DataflowWriteChannel channel(Collection values) {
        final target = CH.emitAndClose(CH.queue(), values)
        NodeMarker.addSourceNode('channel',target)
        return target
    }

    /**
     * INTERNAL ONLY API
     * <p>
     * Add the {@code update} method to an {@code Agent} so that it call implicitly
     * the {@code Agent#updateValue} method
     */
    @CompileDynamic
    static void update(Agent self, Closure message ) {
        assert message != null

        self.send {
            message.call(it)
            updateValue(it)
        }

    }

    static private void checkContext(String method, Object operand) {
        if( operand instanceof ComponentDef && !ExecutionStack.withinWorkflow() )
            throw new IllegalArgumentException("Process invocation are only allowed within a workflow context")
    }

    /**
     * Pipe a channel INTO a process or workflow.
     *
     * @param left
     * @param right
     */
    static Object or(DataflowWriteChannel left, ChainableDef right) {
        checkContext('or', right)
        return right.invoke_o(left)
    }

    /**
     * Pipe a channel INTO an operator.
     *
     * @param left
     * @param right
     */
    static Object or(DataflowWriteChannel left, OpCall right) {
        checkContext('or', right)
        return right.setSource(left).call()
    }

    /**
     * Pipe a channel INTO a closure that defines a custom
     * invocation of a process, workflow, or operator.
     *
     * @param left
     * @param right
     */
    static Object or(DataflowWriteChannel left, Closure right) {
        then(left, right)
    }

    /**
     * Pipe a multi-channel INTO a process or workflow.
     *
     * @param left
     * @param right
     */
    static Object or(ChannelOut left, ChainableDef right) {
        checkContext('or', right)
        return right.invoke_o(left)
    }

    /**
     * Pipe a multi-channel INTO an operator.
     *
     * @param left
     * @param right
     */
    static Object or(ChannelOut left, OpCall right) {
        checkContext('or', right)
        right.setSource(left).call()
    }

    /**
     * Pipe a multi-channel INTO a closure that defines a custom
     * invocation of a process, workflow, or operator.
     *
     * @param left
     * @param right
     */
    static Object or(ChannelOut left, Closure right) {
        then(left, right)
    }

    /**
     * Pipe a process or workflow INTO another process or workflow.
     *
     * @param left
     * @param right
     */
    static Object or(ChainableDef left, ChainableDef right) {
        checkContext('or', left)
        checkContext('or', right)

        def out = left.invoke_a(InvokerHelper.EMPTY_ARGS)

        if( out instanceof DataflowWriteChannel )
            return or((DataflowWriteChannel)out, right)

        if( out instanceof ChannelOut )
            return or((ChannelOut)out, right)

        throw new ScriptRuntimeException("Cannot pipe ${fmtType(out)} with ${fmtType(right)}")
    }

    /**
     * Pipe a process or workflow INTO an operator.
     *
     * @param left
     * @param right
     */
    static Object or(ChainableDef left, OpCall right) {
        checkContext('or', left)
        def out = left.invoke_a(InvokerHelper.EMPTY_ARGS)

        if( out instanceof DataflowWriteChannel )
            return or((DataflowWriteChannel)out, right)

        if( out instanceof ChannelOut )
            return or((ChannelOut)out, right)

        throw new ScriptRuntimeException("Cannot pipe ${fmtType(out)} with ${fmtType(right)}")
    }

    /**
     * Pipe a process or workflow INTO a closure that defines a custom
     * invocation of a process, workflow, or operator.
     *
     * @param left
     * @param right
     */
    static Object or(ChainableDef left, Closure right) {
        checkContext('or', left)
        def out = left.invoke_a(InvokerHelper.EMPTY_ARGS)

        if( out instanceof DataflowWriteChannel )
            return then((DataflowWriteChannel)out, right)

        if( out instanceof ChannelOut )
            return then((ChannelOut)out, right)

        throw new ScriptRuntimeException("Pipeline closure did not return a channel or multi-channel")
    }

    /**
     * Compose a channel with another channel into a multi-channel.
     *
     * @param left
     * @param right
     */
    static ChannelOut and(DataflowWriteChannel left, DataflowWriteChannel right) {
        new ChannelOut(List.of(left, right))
    }

    static ChannelOut and(ChannelOut left, DataflowWriteChannel right) {
        def elements = ChannelOut.spread(left)
        elements << right

        new ChannelOut(elements)
    }

    static ChannelOut and(ChannelOut left, ChannelOut right) {
        def elements = ChannelOut.spread(left)
        elements.addAll(ChannelOut.spread(right))

        new ChannelOut(elements)
    }

    /**
     * Compose a process or workflow with another process or workflow
     * such that they receive the same input channels and their outputs
     * are concatenated.
     *
     * @param left
     * @param right
     */
    static CompositeDef and(ChainableDef left, ChainableDef right) {
        checkContext('and', left)
        checkContext('and', right)
        return new CompositeDef().add(left).add(right)
    }

    static CompositeDef and(CompositeDef left, ChainableDef right) {
        checkContext('and', left)
        checkContext('and', right)
        left.add(right)
    }

}
