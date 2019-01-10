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

package nextflow.script

import nextflow.Channel
import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class EachInParamTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def testNormalize() {

        given:
        def channel = Channel.from(1,2,3,5)
        def value = Channel.value('a')
        def list = Channel.value([4,5,6])
        def each = new EachInParam(Mock(Binding), [])

        expect:
        each.normalizeToVariable(1).val == [1]
        each.normalizeToVariable([3,4,5]).val == [3,4,5]
        each.normalizeToVariable(channel).val == [1,2,3,5]
        each.normalizeToVariable(value).val == ['a']
        each.normalizeToVariable(list).val == [4,5,6]

    }

}
