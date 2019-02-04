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

package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CheckHelperTest extends Specification {

    def 'test is valid' () {

        expect:
        CheckHelper.isValid( 1, Integer )
        CheckHelper.isValid( 1, [1,2,3] )
        CheckHelper.isValid( 1, 1 )
        CheckHelper.isValid( 10, ~/\d+/ )
        !CheckHelper.isValid( 'abc', ~/\d+/ )

        CheckHelper.isValid( 'abc', ~/a*b*c/ )
        !CheckHelper.isValid( 'abz', ~/a*b*c/ )

        !CheckHelper.isValid( 1, [2,3] )
        !CheckHelper.isValid( 1, String )
        !CheckHelper.isValid( 1, 3 )
        !CheckHelper.isValid( 1, null )
        CheckHelper.isValid( null, null )

    }


    def 'test checkParams with map' () {

        when:
        CheckHelper.checkParams('hola', [x:1], [x:[1,2,3]] )
        then:
        true

        when:
        CheckHelper.checkParams('hola', [x:1], [x:[2,3]] )
        then:
        def e = thrown(IllegalArgumentException)
        e.message == "Value '1' cannot be used in in parameter 'x' for operator 'hola' -- Possible values: 2, 3"

        when:
        CheckHelper.checkParams('hola', [x:1], [x:Integer] )
        then:
        true

        when:
        CheckHelper.checkParams('hola', [x:1], [x:String] )
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Value '1' cannot be used in in parameter 'x' for operator 'hola' -- Value don't match: class java.lang.String"

        when:
        CheckHelper.checkParams('hola', [x:true], [x:[Boolean, [1,2,3]]] )
        then:
        true

        when:
        CheckHelper.checkParams('hola', [x:2], [x:[Boolean, [1,2,3]]] )
        then:
        true

        when:
        CheckHelper.checkParams('hola', [x:'Ciao'], [x:[Boolean, [1,2,3]]] )
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Value 'Ciao' cannot be used in in parameter 'x' for operator 'hola' -- Possible values: class java.lang.Boolean, [1, 2, 3]"

    }

    def 'test checkParams splitter map' () {

        def valid = [
                each: Closure,
                by: Integer,
                into: [ Collection ],
                record: [ Boolean, Map ],
                autoClose: Boolean,
                meta: ['file','path','index']
        ]

        when:
        CheckHelper.checkParams ('splitter', [into: [], count: 2], valid)
        then:
        thrown(IllegalArgumentException)

    }


    def 'test checkParams with list' () {

        when:
        CheckHelper.checkParams('hola', [x:1], 'x', 'y' )
        CheckHelper.checkParams('hola', [x:1, y:2], 'x', 'y' )
        then:
        true

        when:
        CheckHelper.checkParams('hola', [z:1], 'x', 'y' )
        then:
        thrown(IllegalArgumentException)

    }

}
