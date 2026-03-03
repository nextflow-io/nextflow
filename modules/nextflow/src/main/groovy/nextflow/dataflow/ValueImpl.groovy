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
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.extension.CH
import nextflow.extension.DataflowHelper

/**
 * Implements the Value type for dataflow v2.
 *
 * @see nextflow.script.types.Value
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ValueImpl {

    private static Session getSession() { Global.getSession() as Session }

    private DataflowVariable source

    ValueImpl(DataflowVariable source) {
        this.source = source
    }

    DataflowVariable getSource() {
        return source
    }

    ChannelImpl flatMap(Function<?,Iterable> transform = null) {
        final target = CH.create()
        final onNext = { value ->
            final iterable = transform != null ? transform.apply(value) : value
            for( final e : iterable )
                target << e
            target << CH.stop()
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("flatMap", [source], [target])
        return new ChannelImpl(target)
    }

    ValueImpl map(Function<?,?> transform) {
        final target = CH.value()
        final onNext = { value ->
            target.bind(transform.apply(value))
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("map", [source], [target])
        return new ValueImpl(target)
    }

    void subscribe(Closure onNext) {
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    void subscribe(Map<String,Closure> events) {
        DataflowHelper.subscribeImpl(source, events)
        NodeMarker.addOperatorNode("subscribe", [source], [])
    }

    ValueImpl view(Function<?,?> transform = null) {
        final target = CH.value()

        final onNext = { value ->
            final result = transform != null ? transform.call(value) : value
            session.printConsole(result?.toString(), true)
            target.bind(value)
        }
        DataflowHelper.subscribeImpl(source, [onNext: onNext])
        NodeMarker.addOperatorNode("view", [source], [target])
        return new ValueImpl(target)
    }

}
