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
