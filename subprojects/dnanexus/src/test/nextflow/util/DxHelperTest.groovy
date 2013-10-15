/*
 * Copyright (c) 2012-2013, the authors.
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
package nextflow.executor
import com.fasterxml.jackson.databind.node.ArrayNode
import nextflow.util.DxHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DxHelperTest extends Specification {


    def testToJsonNode(){

        when:
        def map = [:]
        map.put("one",1); // autoboxed to an object
        map.put("two", [1, 2]); // array of ints is an object
        map.put("three","hello"); // string is an object
        map.put('empty', null )
        map.put('map', [a:10,b:20,c:30])
        map.put('array', ['alpha','beta'] as String[] )
        def result = DxHelper.objToJson(map);

        then:
        result.get('one').intValue() == 1
        result.get('two') instanceof ArrayNode
        result.get('two').get(0).intValue() == 1
        result.get('two').get(1).intValue() == 2
        result.get('three').isTextual()
        result.get('three').textValue() == 'hello'
        result.get('empty').isNull()
        result.get('map').isObject()
        result.get('map').get('a').isInt()
        result.get('map').get('a').intValue() == 10
        result.get('map').get('b').isInt()
        result.get('map').get('b').intValue() == 20
        result.get('map').get('c').isInt()
        result.get('map').get('c').intValue() == 30
        result.get('map').size() == 3
        result.get('array').isArray()
        result.get('array').get(0).isTextual()
        result.get('array').get(0).textValue() == 'alpha'
        result.get('array').get(1).isTextual()
        result.get('array').get(1).textValue() == 'beta'
        result.get('array').size() == 2

    }

    def testJsonNodeToMap2() {

        when:
        def EXPECTED = [field1: "string", field2:[1,2,3], field3:[x: 1, y: 2], field4:[ [p:10, q:20], [w:99], ['a', 'b', 'c'] ]]
        def node = DxHelper.objToJson(EXPECTED)
        def map = DxHelper.jsonToObj(node)

        then:
        map == EXPECTED
        println("MAPA >> $EXPECTED")

    }



    def testJsonNodeToMapNative(){

        when:
        def node = DxHelper.objToJson( [field1: false as Boolean , field2:234.1 as Double, field3: 3 as Integer, field4:2344444444443 as Long, field5:34 as Number, field6:'Hello'])
        def map = DxHelper.jsonToObj(node)

        then:
        map.field1 instanceof Boolean
        map.field2 instanceof Double
        map.field3 instanceof Integer
        map.field4 instanceof Long
        map.field5 instanceof Number
        map.field6 instanceof String
    }

    def testJsonNodeToMapBigDecimal(){

        when:
        def node = DxHelper.objToJson( [field1: 0.23546769425757967900001 as BigDecimal])
        def map = DxHelper.jsonToObj(node)

        then:
        map.field1 == 0.23546769425757967900001 as BigDecimal

        when:
        node = DxHelper.objToJson( [field1: 239423904802384 as BigInteger])
        map = DxHelper.jsonToObj(node)
        then:
        map.field1 == 239423904802384 as BigInteger

        when:
        node = DxHelper.objToJson( [field1: [0xa,0x2,0xf] as byte[] ] )
        map = DxHelper.jsonToObj(node)

        then:
        map.field1 == [0xa,0x2,0xf] as byte[]

    }


    def testJSonNodeToMapArray(){

        when:
        def node = DxHelper.objToJson( [field1:[1,2,3] as Object[] ])
        def map = DxHelper.jsonToObj(node)

        then:
        map.field1 == [1,2,3]
    }


    def testMapToJsonNodeGeneric() {

        when:
        def node = DxHelper.objToJson( [a:1] )
        then:
        node.get('a').intValue() == 1


        when:
        node = DxHelper.objToJson( [a:1, b: 'string'] )
        then:
        node.get('a').intValue() == 1
        node.get('b').textValue() == 'string'


        when:
        node = DxHelper.objToJson( [a:1, b: [1,2,3] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 3
        (node.get('b') as ArrayNode).get(0).intValue() == 1
        (node.get('b') as ArrayNode).get(1).intValue() == 2
        (node.get('b') as ArrayNode).get(2).intValue() == 3


        when:
        node = DxHelper.objToJson( [a:1, b: ['x', 'y'] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 2
        (node.get('b') as ArrayNode).get(0).textValue() == 'x'
        (node.get('b') as ArrayNode).get(1).textValue() == 'y'


        when:
        node = DxHelper.objToJson( [a:1, b: [ [x:1], [y:2, z:3] ] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 2
        (node.get('b') as ArrayNode).get(0).get('x').intValue() == 1


    }


}
