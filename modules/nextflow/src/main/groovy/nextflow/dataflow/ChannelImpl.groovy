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
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.dataflow.ops.CrossOpV2
import nextflow.dataflow.ops.GroupByOp
import nextflow.dataflow.ops.JoinOpV2
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.types.Tuple
import nextflow.util.HashBag

/**
 * Implements the Channel type for dataflow v2.
 *
 * @see nextflow.script.types.Channel
 * @see nextflow.script.types.ChannelV2
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelImpl {

    private static Session getSession() { Global.getSession() as Session }

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

    ChannelImpl cross(Object other) {
        DataflowReadChannel left = this.getReadChannel()
        DataflowReadChannel right
        if( other instanceof ChannelImpl ) {
            right = other.getReadChannel()
        }
        else if( other instanceof ValueImpl ) {
            right = other.getSource()
        }
        else {
            right = CH.value()
            right.bind(other)
        }

        final target = new CrossOpV2(left, right).apply()
        NodeMarker.addOperatorNode("cross", [left, right], [target])
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
        final left = this.getReadChannel()
        final right = other.getReadChannel()
        final target = new JoinOpV2(left, right, opts).apply()
        NodeMarker.addOperatorNode("join", [left, right], [target])
        return new ChannelImpl(target)
    }

    ChannelImpl map(Function<?,?> transform) {
        final source = getReadChannel()
        final target = CH.create()
        final onNext = { value ->
            target << transform.apply(value)
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

    ChannelImpl unique() {
        return unique { v -> v }
    }

    ChannelImpl unique(Function<?,?> transform) {
        final source = getReadChannel()
        final target = CH.create()
        final history = new HashSet<>()

        final onNext = { value ->
            final key = transform.call(value)
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

    ChannelImpl view(Function<?,?> transform = null) {
        final source = getReadChannel()
        final target = CH.create()

        final onNext = { value ->
            final result = transform != null ? transform.call(value) : value
            session.printConsole(result?.toString(), true)
            target << value
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        NodeMarker.addOperatorNode("view", [source], [target])
        return new ChannelImpl(target)
    }

}
