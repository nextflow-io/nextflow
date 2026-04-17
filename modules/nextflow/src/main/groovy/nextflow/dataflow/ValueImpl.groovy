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

import java.util.function.Function

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.dag.NodeMarker
import nextflow.dataflow.ops.CombineOpV2
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.extension.OpCall
import nextflow.extension.OperatorImpl
import nextflow.script.DataflowTypeHelper

/**
 * Implements the Value type for typed workflows.
 *
 * @see nextflow.script.types.Value
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ValueImpl {

    private DataflowVariable source

    ValueImpl(DataflowVariable source) {
        this.source = source
    }

    DataflowVariable getSource() {
        return source
    }

    ValueImpl combine(Object other) {
        DataflowVariable left = this.getSource()

        // determine right-hand source (dataflow value or named args)
        Object right
        if( other instanceof ValueImpl ) {
            right = other.getSource()
        }
        else if( other instanceof Map ) {
            right = other
        }
        else {
            throw new ScriptRuntimeException("Operator `combine` expected a dataflow value or named arguments but received: ${other} [${other.class.simpleName}]")
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

        final target = (DataflowVariable) new CombineOpV2(left, right).apply()
        NodeMarker.addOperatorNode("combine", inputs, [target])
        return new ValueImpl(target)
    }

    ChannelImpl flatMap(Function<?,Iterable> transform = null) {
        final target = CH.create()
        final onNext = { value ->
            final iterable = transform != null ? transform.apply(value) : value
            if( iterable instanceof Tuple )
                throw new ScriptRuntimeException("Operator `flatMap` expected an Iterable but received a tuple: ${iterable}\n")
            for( final e : iterable )
                target << e
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("flatMap", [source], [target])
        return new ChannelImpl(target)
    }

    ValueImpl map(Closure transform) {
        final target = CH.value()
        final onNext = { value ->
            target.bind(transform.call(value))
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("map", [source], [target])
        return new ValueImpl(target)
    }

    ChannelImpl mix(ChannelImpl other) {
        return other.mix(this)
    }

    void subscribe(Closure onNext) {
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    void subscribe(Map<String,Closure> events) {
        DataflowHelper.subscribeImpl(source, events)
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    ValueImpl view(Closure transform = null) {
        return view(Collections.emptyMap(), transform)
    }

    ValueImpl view(Map opts, Closure transform = null) {
        return (ValueImpl) opCall('view', [opts, transform] as Object[])
    }

    // LEGACY OPERATOR FALLBACK

    def methodMissing(String name, Object args) {
        if( name == 'set' || name == 'tap' )
            throw new ScriptRuntimeException("Operator `$name` is not supported in typed workflows -- use an assignment instead")

        if( OperatorImpl.class.getMethods().any { m -> m.name == name } )
            return opCall(name, args)

        throw new MissingMethodException(name, ValueImpl.class, args)
    }

    private Object opCall(String name, Object args) {
        final argsV1 = DataflowTypeHelper.normalizeArray(args, false)
        final result = new OpCall(OperatorImpl.instance, this.source, name, argsV1).call()
        return DataflowTypeHelper.normalizeV2(result)
    }
}
