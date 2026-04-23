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

package nextflow.dataflow

import nextflow.dataflow.ChannelNamespace as channel
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import spock.lang.Specification
import spock.lang.Timeout

import static nextflow.Nextflow.record
import static nextflow.Nextflow.tuple
import static test.ScriptHelper.runDataflow
import static test.ScriptHelper.runScript
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(5)
class ValueImplTest extends Specification {

    def testCombine() {

        when:
        def result = runDataflow {
            channel.value(1).combine(channel.value(2))
        }
        then:
        result.val == tuple(1, 2)
        result.val == tuple(1, 2)  // ValueImpl is re-readable

        when:
        result = runDataflow {
            channel.value(tuple(1, 3)).combine(channel.value(tuple(2, 4)))
        }
        then:
        result.val == tuple(1, 3, 2, 4)

        when:
        runDataflow {
            channel.value(1).combine(channel.of(1, 2, 3))
        }
        then:
        thrown(ScriptRuntimeException)
    }

    def testCombineNamedArgs() {
        when:
        def result = runDataflow {
            channel.value(record(id: 1, a: 'x')).combine(b: channel.value('foo'), c: 'bar')
        }
        then:
        result.val == record(id: 1, a: 'x', b: 'foo', c: 'bar')
    }

    def testFlatMap() {
        when:
        def result = runDataflow {
            channel.value([1,2,3]).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == CH.stop()

        when:
        result = runDataflow {
            channel.value([1,2]).flatMap { v -> [v, v.reverse()] }
        }
        then:
        result.val == [1,2]
        result.val == [2,1]
        result.val == CH.stop()
    }

    def testMap() {
        when:
        def result = runDataflow {
            channel.value('Hello').map { it.reverse() }
        }
        then:
        result.val == 'olleH'
        result.val == 'olleH'
        result.val == 'olleH'
    }

    def testMixWithValue() {
        when:
        def result = runDataflow {
            channel.value(1).mix( channel.of(2,3) ).collect()
        }
        then:
        result.val.sort() == [1,2,3]
    }

    def testSubscribe() {

        when:
        def count = 0
        runDataflow {
            channel.value(42).subscribe { count++ }
        }
        sleep 100
        then:
        count == 1

        when:
        count = 0
        def done = false
        runDataflow {
            channel.value(1).subscribe(
                onNext: { count++ },
                onComplete: { done = true }
            )
        }
        sleep 100
        then:
        done
        count == 1
    }

    def testView() {

        when:
        def result = runDataflow {
            channel.value(42).view()
        }
        then:
        result.val == 42
        result.val == 42  // ValueImpl is re-readable

        when:
        result = runDataflow {
            channel.value(42).view { v -> "value: $v" }
        }
        then:
        result.val == 42
        result.val == 42
    }

    def 'should fall back to legacy operators' () {
        when:
        def result = runScript(
            '''\
            nextflow.enable.types = true

            workflow {
                channel.value([1, 'alpha']).cross( channel.value([1, 'x']) )
            }
            '''
        )
        then:
        result.val == [[1, 'alpha'], [1, 'x']]

        when:
        runDataflow {
            channel.value(1).foo()
        }
        then:
        thrown(MissingMethodException)
    }

}
