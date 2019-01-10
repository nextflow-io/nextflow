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
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RandomSampleTest extends Specification {

    def setupSpec() {
        new Session()
    }

    def 'should produce random sample' () {

        given:
        def ch = Channel.from('A'..'Z')
        def sampler = new RandomSampleOp(ch, 10)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }


    def 'should produce random sample given a short channel' () {

        given:
        def ch = Channel.from('A'..'J')
        def sampler = new RandomSampleOp(ch, 20)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

    def 'should produce random sample given a channel emitting the same number of items as the buffer' () {

        given:
        def ch = Channel.from('A'..'J')
        def sampler = new RandomSampleOp(ch, 10)

        when:
        def result = (List)sampler.apply().toList().val
        then:
        result.size() == 10
        result.unique().size() == 10
        result != 'A'..'J'
    }

    def 'should always produce the same sequence' () {
        given:
        def testSeq = 'A'..'J'
        def ch1 = Channel.from(testSeq)
        def ch2 = Channel.from(testSeq)
        def seed  = 23
        def firstSampler = new RandomSampleOp(ch1, 10, seed)
        def secondSampler = new RandomSampleOp(ch2, 10, seed)

        when:
        def resultFirstRun = (List)firstSampler.apply().toList().val
        def resultSecondRun = (List)secondSampler.apply().toList().val

        then:
        resultFirstRun == resultSecondRun
    }

}
