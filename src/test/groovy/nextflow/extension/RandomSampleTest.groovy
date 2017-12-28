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
