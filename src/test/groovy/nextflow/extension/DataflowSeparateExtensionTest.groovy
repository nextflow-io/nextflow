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
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DataflowSeparateExtensionTest extends Specification{

    @Shared
    Session session

    def setup() {
        session = new Session()
    }

    def cleanup() {
        assert !session.dag.isEmpty()
    }


    def testSeparate() {

        when:
        def str = 'abcdef'
        def (ch1, ch2) = Channel.from(0..3).separate(2) { [it, str[it]] }
        then:
        ch1.val == 0
        ch1.val == 1
        ch1.val == 2
        ch1.val == 3
        ch1.val == Channel.STOP

        ch2.val == 'a'
        ch2.val == 'b'
        ch2.val == 'c'
        ch2.val == 'd'
        ch2.val == Channel.STOP
    }


    def testSeparate2() {

        when:
        def (ch3, ch4) = Channel.from(0..3).map { [it, it+1] } .separate(2)
        then:
        ch3.val == 0
        ch3.val == 1
        ch3.val == 2
        ch3.val == 3
        ch3.val == Channel.STOP

        ch4.val == 1
        ch4.val == 2
        ch4.val == 3
        ch4.val == 4
        ch4.val == Channel.STOP

    }

    def testSeparate3() {

        when:
        def s1 = Channel.create()
        def s2 = Channel.create()
        def s3 = Channel.create()

        Channel.from(1,2,3,4)
                .separate([s1,s2,s3]) { item -> [item+1, item*item, item-1] }

        then:
        s1.val == 2
        s1.val == 3
        s1.val == 4
        s1.val == 5
        s1.val == Channel.STOP
        s2.val == 1
        s2.val == 4
        s2.val == 9
        s2.val == 16
        s2.val == Channel.STOP
        s3.val == 0
        s3.val == 1
        s3.val == 2
        s3.val == 3
        s3.val == Channel.STOP

    }


    def testSeparate4() {
        when:
        def x = Channel.create()
        def y = Channel.create()
        def source = Channel.from([1,2], ['a','b'], ['p','q'])
        source.separate(x,y)
        then:
        x.val == 1
        x.val == 'a'
        x.val == 'p'
        x.val == Channel.STOP
        y.val == 2
        y.val == 'b'
        y.val == 'q'
        y.val == Channel.STOP

        when:
        def x2 = Channel.create()
        def y2 = Channel.create()
        def source2 = Channel.from([1,2], ['a','c','b'], 'z')
        source2.separate(x2,y2)
        then:
        x2.val == 1
        x2.val == 'a'
        x2.val == 'z'
        x2.val == Channel.STOP
        y2.val == 2
        y2.val == 'c'
        y2.val == null
        y2.val == Channel.STOP
    }

}
