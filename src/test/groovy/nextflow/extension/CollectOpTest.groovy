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
import nextflow.util.ArrayBag
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CollectOpTest extends Specification {

    def setupSpec() {
        new Session()
    }

    @Timeout(1)
    def 'should collect items into a list'() {

        when:
        def source = Channel.from(1,2,3)
        def result = source.collect()
        then:
        result.val == [1,2,3]
        result.val instanceof ArrayBag

        when:
        source = Channel.empty()
        result = source.collect()
        then:
        result.val == Channel.STOP

    }

    def 'should collect and invoke a closure on each entry' () {
        given:
        def source = ['hello', 'ciao', 'bonjour']

        when:
        def result = source.channel().collect { it.length() }
        then:
        result.val == [5, 4, 7]
        result.val instanceof ArrayBag
    }

    @Timeout(1)
    def 'should collect and flatten items'() {

        given:
        def source = [[1,['a','b']], [3,['c','d']], [5,['p','q']]]

        when:
        def result = source.channel().collect()
        then:
        result.val == [1,['a','b'],3,['c','d'],5,['p','q']]
        result.val instanceof ArrayBag

        when:
        result = source.channel().collect(flat: true)
        then:
        result.val == [1,['a','b'],3,['c','d'],5,['p','q']]
        result.val instanceof ArrayBag

        when:
        result = source.channel().collect(flat: false)
        then:
        result.val == [[1,['a','b']], [3,['c','d']], [5,['p','q']]]
        result.val instanceof ArrayBag

        when:
        result = source.channel().collect { it.flatten() }
        then:
        result.val == [1,'a','b',3,'c','d',5,'p','q']
        result.val instanceof ArrayBag

    }

    @Timeout(1)
    def 'should collect items into a sorted list '() {

        when:
        def result = [3,1,4,2].channel().collect(sort: true)
        then:
        result.val == [1,2,3,4]
        result.val instanceof ArrayBag

        when:
        result = ['aaa','bb', 'c'].channel().collect(sort: {it->it.size()} as Closure)
        then:
        result.val == ['c','bb','aaa']
        result.val instanceof ArrayBag

        when:
        result = ['aaa','bb', 'c'].channel().collect(sort: {a,b -> a.size()<=>b.size()} as Comparator)
        then:
        result.val == ['c','bb','aaa']
        result.val instanceof ArrayBag

    }


}
