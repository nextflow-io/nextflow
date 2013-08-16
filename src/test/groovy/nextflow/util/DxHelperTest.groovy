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
package nextflow.util

import com.fasterxml.jackson.databind.node.ArrayNode
import com.fasterxml.jackson.databind.node.ObjectNode
import groovy.json.JsonSlurper
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DxHelperTest extends Specification {

    def testDownload() {

        when:
        File target = File.createTempFile('dx-download',null)
        DxHelper.downloadFile( 'file-B7vQg7j0j583xBb7QJV00284', target )

        then:
        target.text.startsWith('>1aboA')

        cleanup:
        target?.delete()

    }


    def testUpload() {

        setup:
        def home = System.getProperty('user.home')
        def config = new File(home, '.dnanexus_config/environment.json')
        def json = new JsonSlurper().parseText(config.text)
        def projectId = json.DX_PROJECT_CONTEXT_ID as String

        when:
        File target = File.createTempFile('dxupload', 'test')
        target.text = "Hello there!\n"
        String id = DxHelper.uploadFile(target, 'test-upload', projectId)

        then:
        id.startsWith('file-')

        cleanup:
        target?.delete()
    }


    def testMapToJsonNode(){

        when:
        Map<String,Object> map = new HashMap<String,Object>();
        map.put("one",1); // autoboxed to an object
        map.put("two", [1, 2]); // array of ints is an object
        map.put("three","hello"); // string is an object
        ObjectNode result = DxHelper.toJsonNode(map);

        then:
        println("RESULTADO >> ${result}")
    }

    def testJsonNodeToMap2() {

        when:
        def node = DxHelper.toJsonNode( [field1: "string", field2:[1,2,3], field3:[x: 1, y: 2], field4:[ [p:10, q:20], [w:99], ['a', 'b', 'c'] ]] )
        def map = DxHelper.jsonNodeToMap(node)

        println map.field2.class.name

        then:
        map.field1 == 'string'
        map.field2 == [1,2,3]
        map.field3 == [x: 1, y: 2]
        map.field4 == [ [p:10, q:20], [w:99], ['a', 'b', 'c'] ]

    }

    def testJsonNodeToMapNative(){

        when:
        def node = DxHelper.toJsonNode( [field1: false as Boolean , field2:234.1 as Double, field3: 3 as Integer, field4:2344444444443 as Long, field5:34 as Number, field6:'Hello'])
        def map = DxHelper.jsonNodeToMap(node)

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
        def node = DxHelper.toJsonNode( [field1: 0.23546769425757967900001 as BigDecimal])
        def map = DxHelper.jsonNodeToMap(node)

        then:
        map.field1 == 0.23546769425757967900001 as BigDecimal
    }

    def testJsonNodeToMapBigInteger(){

        when:
        def node = DxHelper.toJsonNode( [field1: 239423904802384 as BigInteger])
        def map = DxHelper.jsonNodeToMap(node)

        then:
        map.field1 == 239423904802384 as BigInteger
    }

    def testJsonNodeToMapBinary() {

        when:
        def node = DxHelper.toJsonNode( [field1: [0xa,0x2,0xf] as byte[] ] )
        def map = DxHelper.jsonNodeToMap(node)

        then:
        map.field1 == [0xa,0x2,0xf] as byte[]

    }


    def testJSonNodeToMapArray(){

        when:
        def node = DxHelper.toJsonNode( [field1:[1,2,3] as Object[] ])
        def map = DxHelper.jsonNodeToMap(node)

        then:
        map.field1 == [1,2,3]
    }


    def testMapToJsonNodeGeneric() {

        when:
        def node = DxHelper.toJsonNode( [a:1] )
        then:
        node.get('a').intValue() == 1


        when:
        node = DxHelper.toJsonNode( [a:1, b: 'string'] )
        then:
        node.get('a').intValue() == 1
        node.get('b').textValue() == 'string'


        when:
        node = DxHelper.toJsonNode( [a:1, b: [1,2,3] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 3
        (node.get('b') as ArrayNode).get(0).intValue() == 1
        (node.get('b') as ArrayNode).get(1).intValue() == 2
        (node.get('b') as ArrayNode).get(2).intValue() == 3


        when:
        node = DxHelper.toJsonNode( [a:1, b: ['x', 'y'] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 2
        (node.get('b') as ArrayNode).get(0).textValue() == 'x'
        (node.get('b') as ArrayNode).get(1).textValue() == 'y'


        when:
        node = DxHelper.toJsonNode( [a:1, b: [ [x:1], [y:2, z:3] ] ] )
        then:
        node.get('a').intValue() == 1
        node.get('b') instanceof ArrayNode
        (node.get('b') as ArrayNode).size() == 2
        (node.get('b') as ArrayNode).get(0).get('x').intValue() == 1




    }


}
