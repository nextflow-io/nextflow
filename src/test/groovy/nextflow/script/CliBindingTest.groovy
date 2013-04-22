/*
 * Copyright (c) 2012, the authors.
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

package nextflow.script

import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliBindingTest extends Specification {


    def 'test params' () {

        setup:
        def bindings = new CliBinding(new Session([env:[HOME:'/this/path']]))
        bindings.setParams( [field1: 1, field2: 'dos'] )
        bindings.setArgs('a','b','c')

        when:
        // set a generic value
        bindings.setVariable('variable_x', 99)

        // 'params' cannot be overridden
        bindings.setVariable('params', 'xx')
        // 'args' cannot be overridden
        bindings.setVariable('args', 'yy')

        then:
        bindings.getVariable('params').field1 == 1
        bindings.getVariable('params').field2 == 'dos'

        bindings.getVariable('args')  == ['a','b','c']
        bindings.getVariable('variable_x') == 99

        // note: the variable is not defined
        !bindings.hasVariable('HOME')

        // note: BUT it fallback on the local environment
        bindings.getVariable('HOME') == '/this/path'

    }

    def 'test read only map ' () {

        setup:
        CliBinding.ReadOnlyMap map1 = new CliBinding.ReadOnlyMap([:], ['x','y'])
        CliBinding.ReadOnlyMap map2 = new CliBinding.ReadOnlyMap([x:1,y:2,z:3])

        when:
        map1.x = 10
        map1.z = 30

        map2.x = 10
        map2.y = 20
        map2.p = 30
        map2.q = 40


        then:
        map1.x == null
        map1.z == 30

        map2.x == 1
        map2.y == 2
        map2.z == 3
        map2.p == 30
        map2.q == 40


    }

    def 'test hyphenToCamelCase' () {

        expect:
        CliBinding.ParamsMap.hyphenToCamelCase('a') == 'a'
        CliBinding.ParamsMap.hyphenToCamelCase('A') == 'A'
        CliBinding.ParamsMap.hyphenToCamelCase('a-b-c-') == 'aBC'
        CliBinding.ParamsMap.hyphenToCamelCase('aa-bb-cc') == 'aaBbCc'
        CliBinding.ParamsMap.hyphenToCamelCase('alpha-beta-delta') == 'alphaBetaDelta'
        CliBinding.ParamsMap.hyphenToCamelCase('Alpha-Beta-delta') == 'AlphaBetaDelta'

    }

    def 'test camelCaseToHyphen' () {

        expect:
        CliBinding.ParamsMap.camelCaseToHyphen('alphaBetaDelta') == 'alpha-beta-delta'
        CliBinding.ParamsMap.camelCaseToHyphen('AlphaBetaDelta') == 'Alpha-beta-delta'
        CliBinding.ParamsMap.camelCaseToHyphen('Field1') == 'Field1'
        CliBinding.ParamsMap.camelCaseToHyphen('FieldUno') == 'Field-uno'
        CliBinding.ParamsMap.camelCaseToHyphen('FieldUNO') == 'Field-UNO'
        CliBinding.ParamsMap.camelCaseToHyphen('FieldA') == 'Field-A'
        CliBinding.ParamsMap.camelCaseToHyphen('FieldAB') == 'Field-AB'
        CliBinding.ParamsMap.camelCaseToHyphen('FieldAb') == 'Field-ab'

    }

    def 'test put params map' () {

        when:
        def map = new CliBinding.ParamsMap()
        map['alphaBeta'] = 1
        map['alphaBeta'] = 2
        map['alpha-beta'] = 3

        then:
        map['alphaBeta'] == 1
        map['alpha-beta'] == 1


        when:
        map = new CliBinding.ParamsMap()
        map['aaa-bbb-ccc'] = 1
        map['aaaBbbCcc'] = 10
        map['AaaBbbCcc'] = 20

        then:
        map['aaaBbbCcc'] == 1
        map['aaa-bbb-ccc'] == 1


        when:
        map = new CliBinding.ParamsMap()
        map['field1'] = 1
        map['field2'] = 2
        map['Field2'] = 3

        then:
        map['field1'] == 1
        map['field-1'] == null
        map['field2']  == 2
        map['Field2']  == 3

    }



}
