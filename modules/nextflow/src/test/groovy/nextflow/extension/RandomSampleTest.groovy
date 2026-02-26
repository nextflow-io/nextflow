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
import spock.lang.Specification

import static test.ScriptHelper.runDataflow
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RandomSampleTest extends Specification {

    def 'should produce random sample' () {

        when:
        def result = (List) runDataflow {
            Channel.of('A'..'Z').randomSample(10).toList()
        }.val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }


    def 'should produce random sample given a short channel' () {

        when:
        def result = (List) runDataflow {
            Channel.of('A'..'J').randomSample(20).toList()
        }.val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

    def 'should produce random sample given a channel emitting the same number of items as the buffer' () {

        when:
        def result = (List) runDataflow {
            Channel.of('A'..'J').randomSample(10).toList()
        }.val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

    def 'should always produce the same sequence' () {
        given:
        def sequence = 'A'..'J'
        def seed = 23

        when:
        def result1 = (List) runDataflow {
            Channel.of(sequence).randomSample(10, seed).toList()
        }.val
        def result2 = (List) runDataflow {
            Channel.of(sequence).randomSample(10, seed).toList()
        }.val

        then:
        result1 == result2
    }

}
