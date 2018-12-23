/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import nextflow.container.ContainerScriptTokens
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerScriptTokensTest extends Specification {


    def 'should remove single quote and double quotes' () {

        expect:
        ContainerScriptTokens.removeQuotes(null) == null
        ContainerScriptTokens.removeQuotes('Hello') == 'Hello'
        ContainerScriptTokens.removeQuotes('"Hello"') == 'Hello'
        ContainerScriptTokens.removeQuotes("'Hello'") == 'Hello'
        ContainerScriptTokens.removeQuotes("''Hello''") == "'Hello'"
        ContainerScriptTokens.removeQuotes("\"'Hello'\"") == "'Hello'"
        ContainerScriptTokens.removeQuotes("''") == ''
        ContainerScriptTokens.removeQuotes('""') == ''
    }


    def 'should return first string token' () {

        expect:
        ContainerScriptTokens.parse( 'hello' ).image == 'hello'
        ContainerScriptTokens.parse( ' hello ' ).image == 'hello'
        ContainerScriptTokens.parse( 'hello/world --foo --bar\nfoo\nbar' ).image == 'hello/world'
        ContainerScriptTokens.parse( ' \n hello/world\nfoo\nbar ' ).image =='hello/world'

        when:
        ContainerScriptTokens.parse(' ')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should return image, index and variable definitions' () {

        String script

        when:
        script = '''
            #!/bin/bash
            container/name --foo --bar
            '''
        then:
        ContainerScriptTokens.parse( script ).image == 'container/name'
        ContainerScriptTokens.parse( script ).index == 1
        ContainerScriptTokens.parse( script ).variables == [:]
        ContainerScriptTokens.parse( script ).lines == ['#!/bin/bash', 'container/name --foo --bar']

        when:
        script = '''
            #!/bin/bash
            THIS='one'
            THAT="two"

            hello/world --foo --bar
            some
            other
            stuff
            '''
        then:
        ContainerScriptTokens.parse( script ).image ==  'hello/world'
        ContainerScriptTokens.parse( script ).index == 4
        ContainerScriptTokens.parse( script ).variables == [THIS:'one', THAT:'two']
        ContainerScriptTokens.parse( script ).lines == ['#!/bin/bash', "THIS='one'", 'THAT="two"', '', 'hello/world --foo --bar', 'some','other', 'stuff']

    }


}
