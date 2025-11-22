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

package nextflow.script

import spock.lang.Specification
import spock.lang.Unroll

import nextflow.processor.ErrorStrategy
import nextflow.util.Duration
import static nextflow.util.CacheHelper.HashMode
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessConfigTest extends Specification {


    def 'should return defaults' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        expect:
        config.shell ==  ['/bin/bash','-ue']
        config.cacheable
        config.maxRetries == 1
        config.maxErrors == -1
        config.errorStrategy == ErrorStrategy.TERMINATE
    }

    def 'should set properties' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        // setting property using method without brackets
        when:
        config.tag = 'val 1'
        then:
        config.tag == 'val 1'

        // generic value assigned like a 'plain' property
        when:
        config.tag = 99
        then:
        config.tag == 99

    }

    @Unroll
    def 'should set fair directive' () {
        given:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        when:
        config.fair = CONFIG
        then:
        config.getFair() == EXPECTED

        where:
        CONFIG      | EXPECTED
        null        | false
        false       | false
        true        | true
    }

    def 'should parse properties'() {

        when:
        def config = new ProcessConfig( maxDuration:'1h' )
        then:
        config.maxDuration as Duration == Duration.of('1h')
    }


    def 'should not throw MissingPropertyException' () {

        when:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)
        def x = config.hola

        then:
        x == null
        noExceptionThrown()

    }


    def 'should check property existence' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        expect:
        config.containsKey('debug')
        config.containsKey('shell')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null

    }


    def 'should set cache attribute'() {

        when:
        def config = new ProcessConfig(map)
        then:
        config.cacheable == result
        config.isCacheable() == result
        config.getHashMode() == mode

        where:
        result | mode               | map
        true   | HashMode.STANDARD  | [:]
        true   | HashMode.STANDARD  | [cache:true]
        true   | HashMode.STANDARD  | [cache:'yes']
        true   | HashMode.DEEP      | [cache:'deep']
        false  | HashMode.STANDARD  | [cache:false]
        false  | HashMode.STANDARD  | [cache:'false']
        false  | HashMode.STANDARD  | [cache:'off']
        false  | HashMode.STANDARD  | [cache:'no']

    }


    def 'should set ext property' () {

        setup:
        def script = Mock(BaseScript)
        def config = new ProcessConfig(script)

        // setting property using method without brackets
        when:
        config.ext.tool = 'blast'
        config.ext.modules = ['a/1', 'b/2']
        config.ext.command = { "echo $foo" }
        then:
        config.ext.tool == 'blast'
        config.ext.modules ==  ['a/1', 'b/2']

        when:
        def task = config.createTaskConfig()
        task.setContext( [foo: 'Hello'] )
        then:
        task.ext.tool == 'blast'
        task.ext.modules.join(',') == 'a/1,b/2'
        task.ext.command == 'echo Hello'

    }

}
