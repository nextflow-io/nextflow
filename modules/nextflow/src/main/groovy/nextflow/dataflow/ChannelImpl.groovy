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
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.dataflow.ops.CrossOpV2
import nextflow.dataflow.ops.JoinOpV2
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.script.types.Bag
import nextflow.script.types.Tuple
import nextflow.util.ArrayTuple
import nextflow.util.HashBag

/**
 * Implements the Channel type for dataflow v2.
 *
 * @see nextflow.script.types.Channel
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelImpl {

    private static Session getSession() { Global.getSession() as Session }

    private DataflowReadChannel source

    ChannelImpl(source) {
        this.source = read0(source)
    }

    private <T> T read0(source) {
        T result
        if( source instanceof DataflowBroadcast )
            result = (T)CH.getReadChannel(source)
        else if( source instanceof DataflowQueue )
            result = (T)CH.getReadChannel(source)
        else
            result = (T)source
        // track this relationship so that the source channel
        // can be retrieved when rendering the DAG
        if( source != result )
            NodeMarker.addDataflowBroadcastPair(result, source)
        return result
    }

    DataflowReadChannel getSource() {
        return source
    }

    ValueImpl collect() {
        final target = CH.value()
        final result = new HashBag<>()
        final onNext = { value ->
            result.add(value)
        }
        final onComplete = {
            target.bind(result)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return new ValueImpl(target)
    }

    ChannelImpl cross(ChannelImpl right) {
        final target = new CrossOpV2(this.source, right.source).apply()
        return new ChannelImpl(target)
    }

    ChannelImpl filter(Closure condition) {
        final target = CH.create()
        final onNext = { value ->
            if( condition.call(value) )
                target << value
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return new ChannelImpl(target)
    }

    ChannelImpl flatMap(Function<?,Iterable> transform = null) {
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
        return new ChannelImpl(target)
    }

    ChannelImpl groupBy(Function<?,Tuple> transform = null) {
        final target = CH.create()

        final groups = new HashMap<?,Bag<?>>()
        final sizes = new HashMap<?,Integer>()
        final emitted = new HashSet<>()

        final onNext = { value ->
            // obtain key, group size, and value
            final ksv = (transform != null ? transform.apply(value) : value) as List

            def key, size, groupValue
            if( ksv.size() == 2 )
                (key, size, groupValue) = [ ksv[0], -1, ksv[1] ]
            else if( ksv.size() == 3 )
                (key, size, groupValue) = [ ksv[0], ksv[1] as Integer, ksv[2] ]
            else
                throw new ScriptRuntimeException("Operator `groupBy` expected a 3-tuple of (key, size, value) or a 2-tuple of (key, value) but received: ${ksv}")

            if( emitted.contains(key) )
                throw new ScriptRuntimeException("Operator `groupBy` received too many values for grouping key: ${key} (expected ${size})")

            // append value to group
            final group = groups.computeIfAbsent(key, (k) -> new HashBag<>())
            group.add(groupValue)

            // set group size
            if( sizes.containsKey(key) ) {
                if( size != sizes[key] )
                    throw new ScriptRuntimeException("Operator `groupBy` received inconsistent group size for key ${key} -- ${size} != ${sizes[key]}\n")
            }
            else {
                sizes[key] = size
            }

            // emit group when it is complete
            if( group.size() == size ) {
                target << new ArrayTuple<>([key, group])
                groups.remove(key)
                emitted.add(key)
            }
        }
        final onComplete = {
            // emit remaining groups
            groups.each { key, values ->
                if( sizes[key] != -1 )
                    throw new ScriptRuntimeException("Operator `groupBy` received too few values for grouping keys: ${key}")
                target << new ArrayTuple<>([key, values])
            }

            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return new ChannelImpl(target)
    }

    ChannelImpl join(Map opts = [:], ChannelImpl right) {
        final target = new JoinOpV2(this.source, right.source, opts).apply()
        return new ChannelImpl(target)
    }

    ChannelImpl map(Function<?,?> transform) {
        final target = CH.create()
        final onNext = { value ->
            target << transform.apply(value)
        }
        final onComplete = {
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return new ChannelImpl(target)
    }

    ChannelImpl mix(Object... others) {
        final List<DataflowReadChannel> sources = []
        sources.add(this.source)

        for( final other : others ) {
            if( other instanceof ChannelImpl )
                sources.add(other.source)
            else if( other instanceof ValueImpl )
                sources.add(other.source)
            else
                throw new ScriptRuntimeException("Operator `mix` expected a channel or dataflow value but received: ${other} [${other.class.simpleName}]")
        }

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

        return new ChannelImpl(target)
    }

    ValueImpl reduce(BiFunction<?,?,?> accumulator) {
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
        return new ValueImpl(target)
    }

    ValueImpl reduce(Object seed, BiFunction<?,?,?> accumulator) {
        final target = CH.value()

        def result = seed
        final onNext = { value ->
            result = accumulator.apply(result, value)
        }
        final onComplete = {
            target.bind(result)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext, onComplete: onComplete])
        return new ValueImpl(target)
    }

    void subscribe(Closure onNext) {
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
    }

    void subscribe(Map<String,Closure> events) {
        DataflowHelper.subscribeImpl(source, events)
    }

    ChannelImpl unique() {
        return unique { v -> v }
    }

    ChannelImpl unique(Function<?,?> transform) {
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
        return new ChannelImpl(target)
    }

    ChannelImpl until(Closure condition) {
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
        return new ChannelImpl(target)
    }

    ChannelImpl view(Function<?,?> transform = null) {
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
        return new ChannelImpl(target)
    }

}
