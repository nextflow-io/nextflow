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
 */

package nextflow.dataflow

import groovy.transform.CompileStatic
import nextflow.dag.NodeMarker
import nextflow.extension.CH
/**
 * Implements the `channel` namespace, which provides the channel factories
 *
 * @see nextflow.script.namespaces.ChannelNamespace
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ChannelNamespace {

    static ChannelImpl empty() {
        final target = CH.emit(CH.create(), CH.stop())
        NodeMarker.addSourceNode('channel.empty', target)
        return new ChannelImpl(target)
    }

    static ChannelImpl fromList(Collection items) {
        final target = CH.create()
        CH.emitAndClose(target, items as List)
        NodeMarker.addSourceNode('channel.fromList', target)
        return new ChannelImpl(target)
    }

    static ChannelImpl of(Object... items) {
        final target = CH.create()
        final values = new ArrayList<>()
        for( final item : items ) {
            if( item instanceof Range )
                values.addAll(item)
            else
                values.add(item)
        }
        values.add(CH.stop())
        CH.emitValues(target, values)
        NodeMarker.addSourceNode('channel.of', target)
        return new ChannelImpl(target)
    }

    static ChannelImpl topic(String name) {
        final target = CH.topic(name)
        NodeMarker.addSourceNode('channel.topic', target)
        return new ChannelImpl(target)
    }

    static ValueImpl value(Object value) {
        final target = CH.value(value)
        NodeMarker.addSourceNode('channel.value', target)
        return new ValueImpl(target)
    }

}
