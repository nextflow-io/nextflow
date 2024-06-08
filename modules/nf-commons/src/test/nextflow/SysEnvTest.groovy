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

package nextflow

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SysEnvTest extends Specification {

    def 'should get default config' () {
        expect:
        SysEnv.get('HOME')  == System.getenv('HOME')
    }

    def 'should push a new env' () {
        when:
        SysEnv.push([HOME: 'hola'])
        then:
        SysEnv.get('HOME') == 'hola'
        and:
        SysEnv.get().size() == 1

        when:
        SysEnv.push([HOME: 'ciao'])
        then:
        SysEnv.get('HOME') == 'ciao'

        when:
        SysEnv.pop()
        then:
        SysEnv.get('HOME') == 'hola'

        when:
        SysEnv.pop()
        then:
        SysEnv.get('HOME') == System.getenv('HOME')

    }

    @Unroll
    def 'should get bool value' () {
        given:
        SysEnv.push(STATE)

        expect:
        SysEnv.getBool('FOO', DEF) == EXPECTED

        where:
        STATE           | DEF       | EXPECTED
        [:]             | false     | false
        [FOO:'false']   | false     | false
        [FOO:'true']    | false     | true
        and:
        [:]             | true      | true
        [FOO:'false']   | false     | false
        [FOO:'true']    | true      | true

    }
}
