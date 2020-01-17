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

package nextflow.trace


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProgressRecordTest extends Specification {

    def 'should create record' () {
        given:
        def rec = new ProgressRecord(10, 'foo')
        expect:
        with(rec) {
            index == 10
            name == 'foo'
            pending == 0
            submitted == 0
            running == 0
            succeeded == 0
            cached == 0
            failed == 0
            aborted == 0
            stored == 0
            ignored == 0
            retries == 0
            !terminated
            !errored
            loadCpus == 0
            loadMemory == 0
            peakRunning == 0
            peakCpus == 0
            peakMemory == 0
        }
    }

    def 'should get counts' () {
        given:
        def PENDING =1
        def SUBMITTED =2
        def RUNNING =3
        def SUCCEEDED =4
        def FAILED =5
        def CACHED =6
        def STORED =7
        and:
        def rec = new ProgressRecord(10, 'foo')

        when:
        rec.pending =PENDING
        rec.submitted =SUBMITTED
        rec.running =RUNNING
        rec.succeeded =SUCCEEDED
        rec.failed =FAILED
        rec.cached =CACHED
        rec.stored =STORED

        then:
        rec.getCompletedCount() == SUCCEEDED+ FAILED+ CACHED+ STORED
        rec.getTotalCount() == PENDING+ SUBMITTED+ RUNNING + SUCCEEDED+ FAILED+ CACHED+ STORED
    }


}
