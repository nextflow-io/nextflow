/*
 * Copyright 2013-2024, Seqera Labs
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
import nextflow.util.ArrayBag
import spock.lang.Specification
import spock.lang.Timeout

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CollectOpTest extends Specification {

    @Timeout(1)
    def 'should collect items into a list'() {

        when:
        def result = runDataflow {
            Channel.of(1,2,3).collect()
        }
        then:
        result.val == [1,2,3]
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.empty().collect()
        }
        then:
        result.val == Channel.STOP

    }

    def 'should collect and invoke a closure on each entry' () {
        when:
        def result = runDataflow {
            Channel.of('hello', 'ciao', 'bonjour').collect { it.length() }
        }
        then:
        result.val == [5, 4, 7]
        result.val instanceof ArrayBag
    }

    @Timeout(1)
    def 'should collect and flatten items'() {

        given:
        def source = [[1,['a','b']], [3,['c','d']], [5,['p','q']]]

        when:
        def result = runDataflow {
            Channel.fromList(source).collect()
        }
        then:
        result.val == [1,['a','b'],3,['c','d'],5,['p','q']]
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.fromList(source).collect(flat: true)
        }
        then:
        result.val == [1,['a','b'],3,['c','d'],5,['p','q']]
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.fromList(source).collect(flat: false)
        }
        then:
        result.val == [[1,['a','b']], [3,['c','d']], [5,['p','q']]]
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.fromList(source).collect { it.flatten() }
        }
        then:
        result.val == [1,'a','b',3,'c','d',5,'p','q']
        result.val instanceof ArrayBag

    }

    @Timeout(1)
    def 'should collect items into a sorted list '() {

        when:
        def result = runDataflow {
            Channel.of(3,1,4,2).collect(sort: true)
        }
        then:
        result.val == [1,2,3,4]
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.of('aaa','bb', 'c').collect(sort: {v -> v.size()} as Closure)
        }
        then:
        result.val == ['c','bb','aaa']
        result.val instanceof ArrayBag

        when:
        result = runDataflow {
            Channel.of('aaa','bb', 'c').collect(sort: {a,b -> a.size()<=>b.size()} as Comparator)
        }
        then:
        result.val == ['c','bb','aaa']
        result.val instanceof ArrayBag

    }


}
