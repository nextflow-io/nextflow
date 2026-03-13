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

import nextflow.Channel
import nextflow.dataflow.ChannelNamespace as channel
import spock.lang.Specification
import spock.lang.Timeout

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Timeout(10)
class ValueImplTest extends Specification {

    def testFlatMap() {
        when:
        def result = runDataflow {
            channel.value([1,2,3]).flatMap()
        }
        then:
        result.val == 1
        result.val == 2
        result.val == 3
        result.val == Channel.STOP
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

    // def testMixWithValue() {
    //     when:
    //     def result = runDataflow {
    //         channel.value(1).mix( channel.of(2,3) ).collect()
    //     }
    //     then:
    //     result.val.sort() == [1,2,3]
    // }

    def testSubscribe() {
        when:
        def count = 0
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

}
