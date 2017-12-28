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
