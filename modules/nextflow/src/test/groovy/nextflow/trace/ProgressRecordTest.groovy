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


import nextflow.processor.ErrorStrategy
import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import nextflow.util.MemoryUnit
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


    def 'should mark pending' () {
        given:
        def PENDING = 10
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.pending = PENDING
        
        when:
        rec.markPending()
        then:
        rec.pending == PENDING +1
        and:
        rec.submitted == 0
        rec.completedCount == 0
    }

    def 'should mark submitted' () {
        given:
        def PENDING = 10
        def SUBMITTED = 20
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getHashLog() >> HASH }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.pending = PENDING
        rec.submitted = SUBMITTED

        when:
        rec.markSubmitted(task)
        then:
        rec.pending == PENDING -1
        rec.submitted == SUBMITTED +1
        and:
        rec.hash == HASH
    }

    def 'should mark running' () {
        given:
        def RUNNING = 10
        def SUBMITTED = 20
        def CPU = 8
        def MEM = MemoryUnit.of('16 GB')
        def LOAD_CPU = 100
        def LOAD_MEM = 200
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPU
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.submitted = SUBMITTED
        rec.loadCpus = LOAD_CPU
        rec.loadMemory = LOAD_MEM

        when:
        rec.markRunning(task)
        then:
        rec.submitted == SUBMITTED -1
        rec.running == RUNNING +1
        rec.loadCpus == LOAD_CPU + CPU
        rec.loadMemory == LOAD_MEM + MEM.toBytes()
        and:
        rec.peakRunning == rec.running
        rec.peakMemory == rec.loadMemory
        rec.peakCpus == rec.loadCpus
    }

    def 'should mark failed' () {
        given:
        def FAILED = 10
        def RETRIES = 2
        def IGNORED = 3
        def RUNNING = 4
        def SUCCEEDED = 5
        def CPUS = 2
        def MEM = 4.GB
        def LOAD_CPUS = 8
        def LOAD_MEM = 20.GB
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPUS
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes

        when:
        rec.markComplete(task)
        then:
        task.failed >> true
        task.getHashLog() >> HASH
        and:
        rec.failed == FAILED +1
        rec.running == RUNNING -1
        rec.loadCpus == LOAD_CPUS - CPUS
        rec.loadMemory == (LOAD_MEM - MEM).bytes
        and:
        rec.retries == RETRIES
        rec.ignored == IGNORED
        rec.succeeded == SUCCEEDED
        rec.errored
    }

    def 'should mark retried' () {
        given:
        def FAILED = 10
        def RETRIES = 2
        def IGNORED = 3
        def RUNNING = 4
        def SUCCEEDED = 5
        def CPUS = 2
        def MEM = 4.GB
        def LOAD_CPUS = 8
        def LOAD_MEM = 20.GB
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPUS
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes

        when:
        rec.markComplete(task)
        then:
        task.failed >> true
        task.getHashLog() >> HASH
        task.getErrorAction() >> ErrorStrategy.RETRY
        and:
        rec.running == RUNNING -1
        rec.loadCpus == LOAD_CPUS - CPUS
        rec.loadMemory == (LOAD_MEM - MEM).bytes
        and:
        rec.failed == FAILED +1
        rec.retries == RETRIES +1
        and:
        rec.ignored == IGNORED
        rec.succeeded == SUCCEEDED
        !rec.errored
    }

    def 'should mark ignore' () {
        given:
        def FAILED = 10
        def RETRIES = 2
        def IGNORED = 3
        def RUNNING = 4
        def SUCCEEDED = 5
        def CPUS = 2
        def MEM = 4.GB
        def LOAD_CPUS = 8
        def LOAD_MEM = 20.GB
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPUS
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes

        when:
        rec.markComplete(task)
        then:
        task.failed >> true
        task.getHashLog() >> HASH
        task.getErrorAction() >> ErrorStrategy.IGNORE
        and:
        rec.running == RUNNING -1
        rec.loadCpus == LOAD_CPUS - CPUS
        rec.loadMemory == (LOAD_MEM - MEM).bytes
        and:
        rec.failed == FAILED +1
        rec.retries == RETRIES
        rec.ignored == IGNORED +1
        rec.succeeded == SUCCEEDED
        !rec.errored
    }

    def 'should mark aborted' () {
        given:
        def FAILED = 10
        def RETRIES = 2
        def IGNORED = 3
        def RUNNING = 4
        def ABORTED = 5
        def SUCCEEDED = 6
        def CPUS = 2
        def MEM = 4.GB
        def LOAD_CPUS = 8
        def LOAD_MEM = 20.GB
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPUS
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.aborted = ABORTED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes

        when:
        rec.markComplete(task)
        then:
        task.aborted >> true
        task.getHashLog() >> HASH
        and:
        rec.running == RUNNING -1
        rec.loadCpus == LOAD_CPUS - CPUS
        rec.loadMemory == (LOAD_MEM - MEM).bytes
        and:
        rec.aborted == ABORTED +1
        rec.failed == FAILED
        rec.retries == RETRIES
        rec.ignored == IGNORED
        rec.succeeded == SUCCEEDED
    }

    def 'should mark succeeded' () {
        given:
        def FAILED = 10
        def RETRIES = 2
        def IGNORED = 3
        def RUNNING = 4
        def ABORTED = 5
        def SUCCEEDED = 6
        def CPUS = 2
        def MEM = 4.GB
        def LOAD_CPUS = 8
        def LOAD_MEM = 20.GB
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) { getConfig() >> Mock(TaskConfig) {
            getCpus() >> CPUS
            getMemory() >> MEM
        } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.aborted = ABORTED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes

        when:
        rec.markComplete(task)
        then:
        task.aborted >> false
        task.failed >> false
        task.getHashLog() >> HASH
        and:
        rec.running == RUNNING -1
        rec.loadCpus == LOAD_CPUS - CPUS
        rec.loadMemory == (LOAD_MEM - MEM).bytes
        and:
        rec.succeeded == SUCCEEDED +1
        rec.aborted == ABORTED
        rec.failed == FAILED
        rec.retries == RETRIES
        rec.ignored == IGNORED
        !rec.errored
    }

    def 'should mark terminated' () {
        given:
        def rec = new ProgressRecord(10, 'foo')
        and:
        assert !rec.terminated

        when:
        rec.markTerminated()
        then:
        rec.terminated
    }

    def 'should mark cached' () {
        given:
        def CACHED = 10
        def STORED = 20
        and:
        def task = Mock(TaskRun) {
            getHashLog() >> 'XYZ'
        }
        and:
        def trace = Mock(TraceRecord)
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.cached = CACHED
        rec.stored = STORED

        when:
        rec.markCached(task, trace)
        then:
        rec.hash == 'XYZ'
        rec.cached == CACHED +1
        rec.stored == STORED
    }

    def 'should mark store' () {
        given:
        def CACHED = 10
        def STORED = 20
        and:
        def task = Mock(TaskRun)
        and:
        def trace = Mock(TraceRecord)
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.cached = CACHED
        rec.stored = STORED

        when:
        rec.markCached(task, null)
        then:
        rec.cached == CACHED
        rec.stored == STORED +1
        and:
        rec.hash == 'skipped'
    }
}
