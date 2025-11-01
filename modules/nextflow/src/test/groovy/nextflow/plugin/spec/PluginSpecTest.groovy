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

package nextflow.plugin.spec

import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.plugin.extension.Factory
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.Operator
import nextflow.plugin.extension.PluginExtensionPoint
import nextflow.script.dsl.Description
import spock.lang.Specification

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class PluginSpecTest extends Specification {

    def "should generate a plugin spec" () {
        given:
        def extensionPoints = [
            'nextflow.plugin.spec.TestConfig',
            'nextflow.plugin.spec.TestExtension',
        ]
        def spec = new PluginSpec(extensionPoints).build()
        and:
        def definitions = spec.definitions.sort { d -> d.spec.name }

        expect:
        definitions.size() == 4
        and:
        definitions[0] == [
            type: 'ConfigScope',
            spec: [
                name: 'hello',
                description: 'The `hello` scope controls the behavior of the `nf-hello` plugin.',
                children: [
                    [
                        type: 'ConfigOption',
                        spec: [
                            name: 'message',
                            description: 'Message to print to standard output when the plugin is enabled.',
                            type: 'String'
                        ]
                    ]
                ]
            ]
        ]
        definitions[1] == [
            type: 'Factory',
            spec: [
                name: 'helloFactory',
                description: null,
                returnType: 'DataflowWriteChannel',
                parameters: []
            ]
        ]
        definitions[2] == [
            type: 'Operator',
            spec: [
                name: 'helloOperator',
                description: null,
                returnType: 'DataflowWriteChannel',
                parameters: [
                    [
                        name: 'arg0',
                        type: 'DataflowReadChannel'
                    ]
                ]
            ]
        ]
        definitions[3] == [
            type: 'Function',
            spec: [
                name: 'sayHello',
                description: 'Say hello to the given targets',
                returnType: 'void',
                parameters: [
                    [
                        name: 'arg0',
                        type: 'List'
                    ]
                ]
            ]
        ]
    }
}


@ScopeName('hello')
@Description("The `hello` scope controls the behavior of the `nf-hello` plugin.")
class TestConfig implements ConfigScope {

    @ConfigOption
    @Description('''
        Message to print to standard output when the plugin is enabled.
    ''')
    String message
}


class TestExtension extends PluginExtensionPoint {

    @Override
    protected void init(Session session) {
    }

    @Factory
    DataflowWriteChannel helloFactory() {
    }

    @Function
    @Description("Say hello to the given targets")
    void sayHello(List<String> targets) {
    }

    @Operator
    DataflowWriteChannel helloOperator(DataflowReadChannel source) {
    }

}
