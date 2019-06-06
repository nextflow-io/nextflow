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

package nextflow.extension

import static nextflow.util.LoggerHelper.*

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.NF
import nextflow.dag.NodeMarker
import nextflow.exception.ScriptRuntimeException
import nextflow.script.ChainableDef
import nextflow.script.ChannelArrayList
import nextflow.script.CompositeDef
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Implements dataflow channel extension methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelEx {

    /**
     * Assign the {@code source} channel to a global variable with the name specified by the closure.
     * For example:
     * <pre>
     *     Channel.from( ... )
     *            .map { ... }
     *            .set { newChannelName }
     * </pre>
     *
     * @param DataflowReadChannel
     * @param holder A closure that must define a single variable expression
     */
    static void set(DataflowWriteChannel source, Closure holder) {
        final name = CaptureProperties.capture(holder)
        if( !name )
            throw new IllegalArgumentException("Missing name to which set the channel variable")

        if( name.size()>1 )
            throw new IllegalArgumentException("Operation `set` does not allow more than one target name")

        NF.binding.setVariable(name[0], source)
    }

    static void set(ChannelArrayList source, Closure holder) {
        final names = CaptureProperties.capture(holder)
        if( names.size() > source.size() )
            throw new IllegalArgumentException("Operation `set` expects ${names.size()} channels but only ${source.size()} are provided")

        for( int i=0; i<source.size(); i++ ) {
            final ch = source[i]
            final nm = names[i]
            NF.binding.setVariable(nm, ch)
        }
    }

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
     * Creates a channel emitting the entries in the collection to which is applied
     *
     * @param values A list of values to be emitted by the resulting channel
     * @return
     */
    static DataflowWriteChannel channel(Collection values) {
        final target = new DataflowQueue()
        final itr = values.iterator()
        while( itr.hasNext() )
            target.bind(itr.next())
        target.bind(Channel.STOP)
        NodeMarker.addSourceNode('channel',target)
        return target
    }

    /**
     * Close a dataflow queue channel binding a {@link Channel#STOP} item
     *
     * @param source The source dataflow channel to be closed.
     */
    @Deprecated
    static DataflowWriteChannel close(DataflowWriteChannel source) {
        if( NF.isDsl2() )
            throw new DeprecationException("Channel `close` method is not supported any more")
        log.warn "The `close` operator is deprecated -- it will be removed in a future release"
        return ChannelFactory.close0(source)
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
        if( !NF.isDsl2() )
            throw new MissingMethodException(method, operand.getClass())
        //
        //if( !ExecutionStack.withinWorkflow() )
        //    throw new IllegalArgumentException("Process invocation are only allowed within a workflow context")
    }

    /**
     * Implements pipe operation between a channel and a process or a sub-workflow
     *
     * @param left A dataflow channel instance
     * @param right A {@link ChainableDef} object eg. a nextflow process
     * @return The channel resulting the pipe operation
     */
    static Object or(DataflowWriteChannel left, ChainableDef right) {
        checkContext('or', left)
        return right.invoke_o(left)
    }

    /**
     * Implements pipe operation between a channel WITH a operator
     *
     * @param left A {@code DataflowWriteChannel} channel as left operand
     * @param right A {@code OpCall} object representing a operator call as right operand
     * @return The resulting channel object
     */
    static Object or(DataflowWriteChannel left, OpCall right) {
        checkContext('or', left)
        return right.setSource(left).call()
    }

    /**
     * Implements pipe operation between a multi-channels WITH a process or a sub-workflow
     *
     * @param left A {@code ChannelArrayList} multi-channel object as left operand
     * @param right A {@code ChainableDef} object representing a process or sub-workflow call as right operand
     * @return The resulting channel object
     */
    static Object or(ChannelArrayList left, ChainableDef right) {
        checkContext('or', left)
        return right.invoke_o(left)
    }

    /**
     * Implements pipe operation between a multi-channels WITH a operator
     *
     * @param left A {@code ChannelArrayList} multi-channel object as left operand
     * @param right A {@code OpCall} object representing a operator call as right operand
     * @return The resulting channel object
     */
    static Object or(ChannelArrayList left, OpCall right) {
        checkContext('or', left)
        if( right.args.size() )
            throw new ScriptRuntimeException("Process multi-output channel cannot be piped with operator ${right.methodName} for which argument is akready provided")

        right
            .setSource(left[0] as DataflowWriteChannel)
            .setArgs(left[1..-1] as Object[])
            .call()
    }

    /**
     * Implements pipe operation between a process or sub-workflow WITH a operator
     *
     * @param left A {@code ChainableDef} object representing a process or sub-workflow call as left operand
     * @param right A {@code OpCall} object representing a operator call as right operand
     * @return The resulting channel object
     */
    static Object or(ChainableDef left, OpCall right) {
        def out = left.invoke_a(InvokerHelper.EMPTY_ARGS)

        if( out instanceof DataflowWriteChannel )
            return or((DataflowWriteChannel)out, right)

        if( out instanceof ChannelArrayList )
            return or((ChannelArrayList)out, right)

        throw new ScriptRuntimeException("Cannot pipe ${fmtType(out)} with ${fmtType(right)}")
    }

    /**
     * Implements pipe operation between a process or sub-workflow WITH another process or sub-workflow
     *
     * @param left A {@code ChainableDef} object representing a process or sub-workflow call as left operand
     * @param right A {@code ChainableDef} object representing a process or sub-workflow call as right operand
     * @return
     */
    static Object or(ChainableDef left, ChainableDef right) {
        def out = left.invoke_a(InvokerHelper.EMPTY_ARGS)

        if( out instanceof DataflowWriteChannel )
            return or((DataflowWriteChannel)out, right)

        if( out instanceof ChannelArrayList )
            return or((ChannelArrayList)out, right)

        throw new ScriptRuntimeException("Cannot pipe ${fmtType(out)} with ${fmtType(right)}")
    }

    static CompositeDef and(ChainableDef left, ChainableDef right) {
        checkContext('and', left)
        return new CompositeDef().add(left).add(right)
    }

    static CompositeDef and(CompositeDef left, ChainableDef right) {
        checkContext('and', left)
        left.add(right)
    }

}
