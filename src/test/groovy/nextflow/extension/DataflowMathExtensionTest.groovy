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

package nextflow.extension

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowMathExtensionTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should return the min value'() {

        expect:
        Channel.from(4,1,7,5).min().val == 1
        Channel.from("hello","hi","hey").min { it.size() } .val == "hi"
        Channel.from("hello","hi","hey").min { a,b -> a.size()<=>b.size() } .val == "hi"
        Channel.from("hello","hi","hey").min { a,b -> a.size()<=>b.size() } .val == "hi"
        Channel.from("hello","hi","hey").min ({ a,b -> a.size()<=>b.size() } as Comparator) .val == "hi"

    }

    def 'should return the max value'() {
        expect:
        Channel.from(4,1,7,5).max().val == 7
        Channel.from("hello","hi","hey").max { it.size() } .val == "hello"
        Channel.from("hello","hi","hey").max { a,b -> a.size()<=>b.size() } .val == "hello"
        Channel.from("hello","hi","hey").max { a,b -> a.size()<=>b.size() } .val == "hello"
        Channel.from("hello","hi","hey").max ({ a,b -> a.size()<=>b.size() } as Comparator) .val == "hello"

    }

    def 'should return the sum'() {
        expect:
        Channel.from(4,1,7,5).sum().val == 17
        Channel.from(4,1,7,5).sum { it * 2 } .val == 34
        Channel.from( [1,1,1], [0,1,2], [10,20,30] ). sum() .val == [ 11, 22, 33 ]
    }


    def 'should return the mean'() {
        expect:
        Channel.from(10,20,30).mean().val == 20
        Channel.from(10,20,30).mean { it * 2 }.val == 40
        Channel.from( [10,20,30], [10, 10, 10 ], [10, 30, 50]).mean().val == [10, 20, 30]
    }

    def 'should convert string to integers' () {

        expect:
        Channel.value('11').toInteger().val == 11

        when:
        def list = Channel.from('1', '4\n', ' 7 ', '100' )
                .toInteger()
                .toList()
                .getVal()

        then:
        list.size() == 4
        list[0] == 1
        list[1] == 4
        list[2] == 7
        list[3] == 100
        list[0] instanceof Integer
        list[1] instanceof Integer
        list[2] instanceof Integer
        list[3] instanceof Integer
    }


    def 'should convert string to long' () {

        expect:
        Channel.value('33').toLong().val == 33L

        when:
        def list = Channel.from('1', '4\n', ' 7 ', '100' )
                .toLong()
                .toList()
                .getVal()

        then:
        list.size() == 4
        list[0] == 1
        list[1] == 4
        list[2] == 7
        list[3] == 100
        list[0] instanceof Long
        list[1] instanceof Long
        list[2] instanceof Long
        list[3] instanceof Long
    }

    def 'should convert string to float' () {


        expect:
        Channel.value('99.1').toFloat().val == 99.1f

        when:
        def list = Channel.from('1', '4\n', ' 7.5 ', '100.1' )
                .toFloat()
                .toList()
                .getVal()

        then:
        list.size() == 4
        list[0] == 1
        list[1] == 4
        list[2] == 7.5f
        list[3] == 100.1f
        list[0] instanceof Float
        list[1] instanceof Float
        list[2] instanceof Float
        list[3] instanceof Float
    }

    def 'should convert string to double' () {

        expect:
        Channel.value('99.1').toDouble().val == 99.1d

        when:
        def list = Channel.from('1', '4\n', ' 7.5 ', '100.1' )
                .toDouble()
                .toList()
                .getVal()

        then:
        list.size() == 4
        list[0] == 1
        list[1] == 4
        list[2] == 7.5d
        list[3] == 100.1d
        list[0] instanceof Double
        list[1] instanceof Double
        list[2] instanceof Double
        list[3] instanceof Double

    }

    def 'should return a random sample' () {

        when:
        def result = Channel
                .from(0,1,2,3,4,5,6,7,8,9)
                .randomSample(5)
                .toList().val as List

        then:
        result.size() == 5
        result.unique().size() == 5
        result != [0,1,2,3,4]
        result[0] in 0..9
        result[1] in 0..9
        result[2] in 0..9
        result[3] in 0..9
        result[4] in 0..9

    }

}
