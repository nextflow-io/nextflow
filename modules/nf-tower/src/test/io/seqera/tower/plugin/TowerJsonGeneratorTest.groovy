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

package io.seqera.tower.plugin

import groovy.json.JsonGenerator
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerJsonGeneratorTest extends Specification {

    def 'should chomp too long values' () {
        given:
        def scheme = [foo: 5, 'bar.one': 5]
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), scheme)

        when:
        def x = gen.toJson( [foo: "Hola", bar: 'mundo'] )
        then:
        x == '{"foo":"Hola","bar":"mundo"}'

        when:
        x = gen.toJson( [foo: "Hello world"] )
        then:
        x == '{"foo":"Hello"}'

        when:
        x = gen.toJson( [bar: [one: "Hello world", two: "Hola mundo"]] )
        then:
        x == '{"bar":{"one":"Hello","two":"Hola mundo"}}'
    }

    def 'should normalise gitmodules attribute' () {
        given:
        def scheme = ['workflow.manifest.gitmodules': 10]
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), scheme)

        when:
        def json = gen.toJson( [workflow: [manifest: [gitmodules: ['a','b','c']]]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"a,b,c"}}}'

        when:
        json = gen.toJson( [workflow: [manifest: [gitmodules: 'abc']]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"abc"}}}'

        when:
        json = gen.toJson( [workflow: [manifest: [gitmodules: '123456789012345']]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"1234567890"}}}'
    }

}
