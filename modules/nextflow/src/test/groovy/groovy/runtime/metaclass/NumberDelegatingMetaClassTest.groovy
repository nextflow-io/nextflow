package groovy.runtime.metaclass
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

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NumberDelegatingMetaClassTest extends Specification {


    def 'should parse duration units' () {
        expect:
        1.ms  == Duration.of('1 millis')
        1.sec == Duration.of('1 sec')
        2.min  == Duration.of('2 min')
        3.hour == Duration.of('3 hours')
        4.days == Duration.of('4 days')
    }

    def 'should multiply number with a duration unit' () {

        expect:
        2.sec * 5 == Duration.of('10 sec')
        5 * 3.sec == Duration.of('15 sec')

    }

    def 'should multiply number with a memory unit'() {
        expect:
        100.KB * 2 == MemoryUnit.of('200 KB')
        10 * 100.KB == MemoryUnit.of('1000 KB')
    }


    def 'should parse memory units' () {
        expect:
        1.kb  == MemoryUnit.of('1 KB')
        2.mb == MemoryUnit.of('2 MB')
        3.gb  == MemoryUnit.of('3 GB')
        4.tb == MemoryUnit.of('4 TB')
        5.pb == MemoryUnit.of('5 PB')

        3.5.GB == MemoryUnit.of('3.5 GB')
    }


}
