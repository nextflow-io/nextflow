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

package nextflow.script

import nextflow.util.ReadOnlyMap
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptBindingTest extends Specification {


    def 'test params' () {

        setup:
        def bindings = new ScriptBinding(env: [HOME:'/this/path'])
        bindings.setParams( [field1: 1, field2: 'dos'] )
        bindings.setArgs(['a','b','c'])

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

        // note: BUT it fallback on the local environment
        bindings.getVariable('HOME') == '/this/path'

    }

    def 'test read only map' () {

        setup:
        ReadOnlyMap map1 = new ReadOnlyMap([:], ['x','y'])
        ReadOnlyMap map2 = new ReadOnlyMap([x:1,y:2,z:3])

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

    def 'should convert hyphen separated string to camel case' () {

        expect:
        ScriptBinding.ParamsMap.hyphenToCamelCase('a') == 'a'
        ScriptBinding.ParamsMap.hyphenToCamelCase('A') == 'A'
        ScriptBinding.ParamsMap.hyphenToCamelCase('a-b-c-') == 'aBC'
        ScriptBinding.ParamsMap.hyphenToCamelCase('aa-bb-cc') == 'aaBbCc'
        ScriptBinding.ParamsMap.hyphenToCamelCase('alpha-beta-delta') == 'alphaBetaDelta'
        ScriptBinding.ParamsMap.hyphenToCamelCase('Alpha-Beta-delta') == 'AlphaBetaDelta'

    }

    def 'should convert camel case string to hyphen separated' () {

        expect:
        ScriptBinding.ParamsMap.camelCaseToHyphen('alphaBetaDelta') == 'alpha-beta-delta'
        ScriptBinding.ParamsMap.camelCaseToHyphen('AlphaBetaDelta') == 'Alpha-beta-delta'
        ScriptBinding.ParamsMap.camelCaseToHyphen('Field1') == 'Field1'
        ScriptBinding.ParamsMap.camelCaseToHyphen('FieldUno') == 'Field-uno'
        ScriptBinding.ParamsMap.camelCaseToHyphen('FieldUNO') == 'Field-UNO'
        ScriptBinding.ParamsMap.camelCaseToHyphen('FieldA') == 'Field-A'
        ScriptBinding.ParamsMap.camelCaseToHyphen('FieldAB') == 'Field-AB'
        ScriptBinding.ParamsMap.camelCaseToHyphen('FieldAb') == 'Field-ab'

    }

    def 'should put an entry in the params map' () {

        when:
        def map = new ScriptBinding.ParamsMap()
        map['alphaBeta'] = 1
        map['alphaBeta'] = 2
        map['alpha-beta'] = 3

        then:
        map['alphaBeta'] == 1
        map['alpha-beta'] == 1


        when:
        map = new ScriptBinding.ParamsMap()
        map['aaa-bbb-ccc'] = 1
        map['aaaBbbCcc'] = 10
        map['AaaBbbCcc'] = 20

        then:
        map['aaaBbbCcc'] == 1
        map['aaa-bbb-ccc'] == 1


        when:
        map = new ScriptBinding.ParamsMap()
        map['field1'] = 1
        map['field2'] = 2
        map['Field2'] = 3

        then:
        map['field1'] == 1
        map['field-1'] == null
        map['field2']  == 2
        map['Field2']  == 3

    }

    def 'should wrap a string value with quote character' () {

        expect:
        ScriptBinding.ParamsMap.wrap('hello',null) == 'hello'
        ScriptBinding.ParamsMap.wrap('hello','"') == '"hello"'
        ScriptBinding.ParamsMap.wrap('hello world',null) == '"hello world"'
        ScriptBinding.ParamsMap.wrap('hello world',"'") == "'hello world'"
        ScriptBinding.ParamsMap.wrap('hello"world',"'") == "'hello\"world'"
        ScriptBinding.ParamsMap.wrap('hello"world',null) == '"hello\\"world"'

    }

    def 'should return a command line formatted string'() {

        when:
        def params = new ScriptBinding.ParamsMap('foo-bar':1)
        then:
        params.size() == 2
        params.fooBar == 1
        params.'foo-bar' == 1
        params.all() == '--foo-bar 1'

        expect:
        new ScriptBinding.ParamsMap(x:1).all() == '--x 1'
        new ScriptBinding.ParamsMap(x:1, y: 2).all() == '--x 1 --y 2'
        new ScriptBinding.ParamsMap(x:1, y: 2).all(sep:'=') == '--x=1 --y=2'
        new ScriptBinding.ParamsMap(x:1, y: 2).all(sep:'=', prefix:'-') == '-x=1 -y=2'
        new ScriptBinding.ParamsMap(x:1, y: 'hello world').all() == '--x 1 --y "hello world"'
        new ScriptBinding.ParamsMap(x:1, y: 'hello world').all(quote:"'") == '--x \'1\' --y \'hello world\''
        new ScriptBinding.ParamsMap(x:1, y: "O'Connors").all(quote:"'") == "--x '1' --y 'O\\'Connors'"

    }

    def 'should get the variable name giving the value'() {

        given:
        def X=1
        def Y='hello'
        def binding = new ScriptBinding([:])
        binding.setVariable('foo', X)
        binding.setVariable('bar', Y)

        expect:
        binding.getVariableName(X) == 'foo'
        binding.getVariableName(Y) == 'bar'
        binding.getVariableName(new String(Y)) == null
        binding.getVariableName(null) == null
    }

}
