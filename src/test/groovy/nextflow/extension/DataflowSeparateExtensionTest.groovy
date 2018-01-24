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
