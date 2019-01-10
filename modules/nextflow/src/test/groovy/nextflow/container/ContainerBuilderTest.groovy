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

package nextflow.container

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerBuilderTest extends Specification {

    def 'should return mount flags'() {

        given:
        def builder = Spy(ContainerBuilder)

        expect:
        builder.mountFlags(false) == ''
        builder.mountFlags(true) == ':ro'

    }

    def 'should make env var' () {
        given:
        StringBuilder result
        def builder = Spy(ContainerBuilder)

        when:
        result = builder.makeEnv([FOO: 'x', BAR: 'y'])
        then:
        result.toString() == '-e "FOO=x" -e "BAR=y"'

        when:
        result = builder.makeEnv('FOO=hello')
        then:
        result.toString() == '-e "FOO=hello"'
        
        when:
        result = builder.makeEnv( 'FOO' )
        then:
        result.toString() == '${FOO:+-e "FOO=$FOO"}'

        when:
        builder.makeEnv( 1 )
        then:
        thrown(IllegalArgumentException)

    }

}
