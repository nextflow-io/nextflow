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

package nextflow.extension.plugin

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.plugin.extension.PluginExtensionProvider
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginExtensionProviderTest extends Specification {

    def 'should check dataflow read channel' () {
        expect:
        PluginExtensionProvider.isReadChannel(DataflowVariable.class)
        PluginExtensionProvider.isReadChannel(DataflowQueue.class)
        !PluginExtensionProvider.isReadChannel(DataflowBroadcast.class)
    }

    def 'should check extension method' () {
        given:
        def ext = PluginExtensionProvider.INSTANCE()
        expect:
        ext.isExtensionMethod(new DataflowVariable(), 'map')
        ext.isExtensionMethod(new DataflowVariable(), 'flatMap')
        !ext.isExtensionMethod(new DataflowVariable(), 'foo')
    }

    def 'should invoke ext method' () {
        given:
        def ext = PluginExtensionProvider.INSTANCE()
        def ch = new DataflowQueue(); ch<<1<<2<<3

        when:
        def result = ext.invokeExtensionMethod(ch, 'map', { it -> it * it })
        then:
        result instanceof DataflowReadChannel
        result.val == 1
        result.val == 4
        result.val == 9
    }
}
