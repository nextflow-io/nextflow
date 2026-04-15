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

package nextflow.dataflow

import java.util.concurrent.atomic.AtomicInteger
import java.util.function.BiFunction
import java.util.function.Function

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.dag.NodeMarker
import nextflow.dataflow.ops.CombineOpV2
import nextflow.dataflow.ops.GroupByOp
import nextflow.dataflow.ops.JoinOpV2
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.extension.OpCall
import nextflow.extension.OperatorImpl
import nextflow.script.DataflowTypeHelper
import nextflow.script.types.Tuple
import nextflow.util.HashBag
import nextflow.util.RecordMap

/**
 * Implements the Channel type for typed workflows.
 *
 * @see nextflow.script.types.Channel
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelImpl {

    private DataflowWriteChannel source

    ChannelImpl(DataflowWriteChannel source) {
        this.source = source
    }

    DataflowWriteChannel getSource() {
        return source
    }

    DataflowReadChannel getReadChannel() {
        return read0(source)
    }

    private static DataflowReadChannel read0(DataflowWriteChannel source) {
        DataflowReadChannel result
        if( source instanceof DataflowBroadcast )
            result = CH.getReadChannel(source)
        else
            throw new IllegalArgumentException()
        // track this relationship so that the source channel
        // can be retrieved when rendering the DAG
        if( source != result )
            NodeMarker.addDataflowBroadcastPair(result, source)
        return result
    }

    ValueImpl collect() {
        final source = getReadChannel()
        final target = CH.value()
        final result = new HashBag<>()
        final onNext = { value ->
            result.add(value)
        }
        final onComplete = {
            target.bind(result)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("collect", [source], [target])
        return new ValueImpl(target)
    }

    ChannelImpl combine(Object other) {
        DataflowReadChannel left = this.getReadChannel()

        // determine right-hand source (channel, dataflow value, or named args)
        Object right
        if( other instanceof ChannelImpl ) {
            right = other.getReadChannel()
        }
        else if( other instanceof ValueImpl ) {
            right = other.getSource()
        }
        else if( other instanceof Map ) {
            right = other
        }
        else {
            throw new ScriptRuntimeException("Operator `combine` expected a channel, dataflow value, or named arguments, but received: ${other} [${other.class.simpleName}]")
        }

        // record dataflow inputs
        final inputs = []
        inputs.add(left)
        if( right instanceof Map ) {
            for( final value : right.values() ) {
                if( value instanceof ValueImpl )
                    inputs.add(value.getSource())
            }
        }
        else {
            inputs.add(right)
        }

        final target = new CombineOpV2(left, right).apply()
        NodeMarker.addOperatorNode("combine", inputs, [target])
        return new ChannelImpl(target)
    }

    ChannelImpl filter(Closure condition) {
        final source = getReadChannel()
        final target = CH.create()
        final onNext = { value ->
            if( condition.call(value) )
                target << value
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("filter", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl flatMap(Function<?,Iterable> transform = null) {
        final source = getReadChannel()
        final target = CH.create()
        final onNext = { value ->
            final iterable = transform != null ? transform.apply(value) : value
            if( iterable instanceof Tuple )
                throw new ScriptRuntimeException("Operator `flatMap` expected an Iterable but received a tuple: ${iterable}\n")
            for( final e : iterable )
                target << e
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("flatMap", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl groupBy() {
        final source = getReadChannel()
        final target = new GroupByOp(source).apply()
        NodeMarker.addOperatorNode("groupBy", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl join(Map opts = [:], ChannelImpl other) {
        if( opts.by instanceof Integer )
            return (ChannelImpl) opCall('join', [opts, other] as Object[])
        final left = this.getReadChannel()
        final right = other.getReadChannel()
        final target = new JoinOpV2(left, right, opts).apply()
        NodeMarker.addOperatorNode("join", [left, right], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl map(Closure transform) {
        final source = this.getReadChannel()
        final target = CH.create()
        final onNext = { value ->
            target << transform.call(value)
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("map", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl mix(Object other) {
        DataflowReadChannel left = this.getReadChannel()
        DataflowReadChannel right
        if( other instanceof ChannelImpl )
            right = other.getReadChannel()
        else if( other instanceof ValueImpl )
            right = other.getSource()
        else
            throw new ScriptRuntimeException("Operator `mix` expected a channel or dataflow value but received: ${other} [${other.class.simpleName}]")

        final sources = [left, right]
        final target = CH.create()
        final count = new AtomicInteger(sources.size())
        final onNext = { value ->
            target << value
        }
        final onComplete = {
            if( count.decrementAndGet() == 0 )
                target << CH.stop()
        }
        for( final ch : sources )
            DataflowHelper.subscribeImpl(ch, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("mix", sources, [target])
        return new ChannelImpl(target)
    }

    ValueImpl reduce(BiFunction<?,?,?> accumulator) {
        final source = getReadChannel()
        final target = CH.value()

        boolean first = true
        def result
        final onNext = { value ->
            if( first ) {
                result = value
                first = false
            }
            else {
                result = accumulator.apply(result, value)
            }
        }
        final onComplete = {
            if( first ) {
                def e = new ScriptRuntimeException("Operator `reduce` received an empty channel with no initial value -- make sure to provide an initial value if the channel might be empty")
                target.bindError(e)
                throw e
            }
            target.bind(result)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("reduce", [source], [target])
        return new ValueImpl(target)
    }

    ValueImpl reduce(Object seed, BiFunction<?,?,?> accumulator) {
        final source = getReadChannel()
        final target = CH.value()

        def result = seed
        final onNext = { value ->
            result = accumulator.apply(result, value)
        }
        final onComplete = {
            target.bind(result)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("reduce", [source], [target])
        return new ValueImpl(target)
    }

    void subscribe(Closure onNext) {
        final source = getReadChannel()
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    void subscribe(Map<String,Closure> events) {
        final source = getReadChannel()
        DataflowHelper.subscribeImpl(source, events)
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    ChannelImpl unique(Closure transform = null) {
        final source = getReadChannel()
        final target = CH.create()
        final history = new HashSet<>()

        final onNext = { value ->
            final key = transform != null ? transform.call(value) : value
            if( !history.contains(key) ) {
                history.add(key)
                target << value
            }
        }
        final onComplete = {
            history.clear()
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("unique", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl until(Closure condition) {
        final source = getReadChannel()
        final target = CH.create()

        boolean done = false
        final onNext = { value ->
            if( done )
                return
            if( condition.call(value) ) {
                target << CH.stop()
                done = true
            }
            else {
                target << value
            }
        }
        final onComplete = {
            if( !done )
                target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("until", [source], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl view(Closure transform = null) {
        return view(Collections.emptyMap(), transform)
    }

    ChannelImpl view(Map opts, Closure transform = null) {
        return (ChannelImpl) opCall('view', [opts, transform] as Object[])
    }

    // LEGACY OPERATOR FALLBACK

    def methodMissing(String name, Object args) {
        if( name == 'set' || name == 'tap' )
            throw new ScriptRuntimeException("Operator `$name` is not supported in typed workflows -- use an assignment instead")

        if( OperatorImpl.class.getMethods().any { m -> m.name == name } )
            return opCall(name, args)

        throw new MissingMethodException(name, ChannelImpl.class, args)
    }

    private Object opCall(String name, Object args) {
        final argsV1 = DataflowTypeHelper.normalizeArray(args, false)
        final result = new OpCall(OperatorImpl.instance, this.source, name, argsV1).call()
        return DataflowTypeHelper.normalizeV2(result)
    }
}
