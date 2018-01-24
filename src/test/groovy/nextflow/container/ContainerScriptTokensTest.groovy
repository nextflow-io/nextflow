/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
