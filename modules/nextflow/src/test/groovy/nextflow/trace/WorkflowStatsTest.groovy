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
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowStatsTest extends Specification {


    def 'should return compute time string' () {

        given:
        WorkflowStats stats

        when:
        stats = new WorkflowStats(succeedMillis: 3_600_000)
        then:
        stats.getComputeTimeFmt() == '1.0'

        when:
        stats = new WorkflowStats(succeedMillis: 10_000_000)
        then:
        stats.getComputeTimeFmt() == '2.8'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, cachedMillis: 20_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (20% cached)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 20_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (20% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 80_000_000, failedMillis: 5_000_000, cachedMillis: 15_000_000)
        then:
        stats.getComputeTimeFmt() == '27.8 (15% cached, 5% failed)'

        when:
        stats = new WorkflowStats(succeedMillis: 180_000)
        then:
        stats.getComputeTimeFmt() == '0.1'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000)
        then:
        stats.getComputeTimeFmt() == '(a few seconds)'

        when:
        stats = new WorkflowStats(succeedMillis: 120_000_000_000)
        then:
        stats.getComputeTimeFmt() == "33'333.3"
    }

    def 'should return cpu time' () {
        given:
        long seconds
        def stats = new WorkflowStats()

        when:
        seconds = stats.getCpuTime( new TraceRecord() )
        then:
        seconds == 0

        when:
        seconds = stats.getCpuTime( new TraceRecord([realtime: 35_000, cpus: 1]) )
        then:
        seconds == 35_000

        when:
        seconds = stats.getCpuTime( new TraceRecord([realtime: 10_000, cpus: 4]) )
        then:
        seconds == 40_000
    }


    def 'should get cpu time' () {
        given:
        def stats = new WorkflowStats()
        def record = Mock(TraceRecord)
        long result

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 1_000
        1 * record.get('cpus') >> 1
        result == 1000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 2_000
        1 * record.get('cpus') >> 4
        result == 8_000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> 2_000
        1 * record.get('cpus') >> null
        result == 2_000

        when:
        result = stats.getCpuTime(record)
        then:
        1 * record.get('realtime') >> null
        1 * record.get('cpus') >> null
        result == 0
    }

    def 'should return formatted counts' () {

        given:
        def stats = new WorkflowStats()

        when:
        stats.succeededCount = 12345
        then:
        stats.succeedCount == 12345
        stats.succeedCountFmt == "12'345"

        when:
        stats.cachedCount = 88774411
        then:
        stats.cachedCount == 88774411
        stats.cachedCountFmt == "88'774'411"

        when:
        stats.ignoredCount = 66332211
        then:
        stats.ignoredCount == 66332211
        stats.ignoredCountFmt == "66'332'211"

        when:
        stats.failedCount = 33776644
        then:
        stats.failedCount == 33776644
        stats.failedCountFmt == "33'776'644"
    }


    def 'should return task percents' () {
        given:
        def stats = new WorkflowStats(succeededCount: 20, cachedCount: 40, ignoredCount: 60, failedCount: 80)

        expect: 
        stats.getSucceedPct() == 10.0f
        stats.getCachedPct() == 20.0f
        stats.getIgnoredPct() == 30.0f
        stats.getFailedPct() == 40.0f

    }

    def 'should clone stats'  () {
        given:
        def stats = new WorkflowStats()
        when:
        def copy = stats.clone()
        then:
        stats == copy

        when:
        stats.markCreated( Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }  )
        and:
        copy = stats.clone()
        then:
        copy == stats
        and:
        copy.progressLength == 1
        copy.processes.size() == 1
        copy.processes[0] == stats.processes[0]
        !copy.processes[0] .is (stats.processes[0])
    }

    def 'should mark pending' () {
        given:
        def PENDING = 10
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.pending = PENDING
        and:
        def stats = new WorkflowStats(
                        pendingCount: PENDING,
                        records: [0: rec])
        
        when:
        stats.markPending( Mock(TaskProcessor) { getId() >> 0 } )
        then:
        stats.pendingCount == PENDING +1
        and:
        rec.pending == PENDING +1
    }

    def 'should mark submitted' () {
        given:
        def PENDING = 10
        def SUBMITTED = 20
        def HASH = 'xyz'
        and:
        def task = Mock(TaskRun) {
            getHashLog() >> HASH
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
        }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.pending = PENDING
        rec.submitted = SUBMITTED
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                pendingCount: PENDING,
                submittedCount: SUBMITTED )

        when:
        stats.markSubmitted(task)
        then:
        stats.pendingCount == PENDING -1
        stats.submittedCount == SUBMITTED +1
        and:
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
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPU; getMemory() >> MEM }
        }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.running = RUNNING
        rec.submitted = SUBMITTED
        rec.loadCpus = LOAD_CPU
        rec.loadMemory = LOAD_MEM
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                submittedCount: SUBMITTED,
                loadCpus: LOAD_CPU,
                loadMemory: LOAD_MEM)

        when:
        stats.markRunning(task)
        then:
        stats.submittedCount == SUBMITTED -1
        stats.runningCount == RUNNING +1
        stats.loadCpus == LOAD_CPU + CPU
        stats.loadMemory == LOAD_MEM + MEM.toBytes()
        and:
        stats.peakRunning == stats.runningCount
        stats.peakMemory == stats.loadMemory
        stats.peakCpus == stats.loadCpus

        and:
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
        def DURATION = Duration.of('5 sec')
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPUS; getMemory() >> MEM }
        }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                failedCount: FAILED,
                retriesCount: RETRIES,
                ignoredCount: IGNORED,
                submittedCount: SUCCEEDED,
                loadCpus: LOAD_CPUS,
                loadMemory: LOAD_MEM.bytes)
        and:
        def trace = Mock(TraceRecord)

        when:
        stats.markCompleted(task, trace)

        then:
        1 * trace.get('realtime')  >> DURATION.millis
        1 * trace.get('cpus') >> 1
        
        and:
        stats.failedCount == FAILED +1
        stats.runningCount == RUNNING -1
        stats.loadCpus == LOAD_CPUS - CPUS
        stats.loadMemory == (LOAD_MEM - MEM).bytes
        stats.failedDuration == DURATION
        and:
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
        def DURATION = 3.sec
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPUS; getMemory() >> MEM } }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                failedCount: FAILED,
                retriesCount: RETRIES,
                ignoredCount: IGNORED,
                submittedCount: SUCCEEDED,
                loadCpus: LOAD_CPUS,
                loadMemory: LOAD_MEM.bytes)
        and:
        def trace = Mock(TraceRecord)

        when:
        stats.markCompleted(task, trace)
        then:
        task.failed >> true
        task.getHashLog() >> HASH
        task.getErrorAction() >> ErrorStrategy.RETRY
        1 * trace.get('realtime')  >> DURATION.millis
        1 * trace.get('cpus') >> 1

        and:
        stats.runningCount == RUNNING -1
        stats.loadCpus == LOAD_CPUS - CPUS
        stats.loadMemory == (LOAD_MEM - MEM).bytes
        stats.failedCount == FAILED +1
        stats.retriesCount == RETRIES +1
        stats.failedDuration == DURATION
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
        def DURATION = 3.sec
        and:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPUS; getMemory() >> MEM } }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                failedCount: FAILED,
                retriesCount: RETRIES,
                ignoredCount: IGNORED,
                succeededCount: SUCCEEDED,
                loadCpus: LOAD_CPUS,
                loadMemory: LOAD_MEM.bytes)
        and:
        def trace = Mock(TraceRecord)

        when:
        stats.markCompleted(task, trace)
        then:
        task.failed >> true
        task.getHashLog() >> HASH
        task.getErrorAction() >> ErrorStrategy.IGNORE
        1 * trace.get('realtime')  >> DURATION.millis
        1 * trace.get('cpus') >> 1
        and:
        stats.runningCount == RUNNING -1
        stats.loadCpus == LOAD_CPUS - CPUS
        stats.loadMemory == (LOAD_MEM - MEM).bytes
        stats.failedCount == FAILED +1
        stats.retriesCount == RETRIES
        stats.ignoredCount == IGNORED +1
        stats.succeededCount == SUCCEEDED
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
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPUS; getMemory() >> MEM } }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.aborted = ABORTED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                failedCount: FAILED,
                retriesCount: RETRIES,
                ignoredCount: IGNORED,
                succeededCount: SUCCEEDED,
                abortedCount: ABORTED,
                loadCpus: LOAD_CPUS,
                loadMemory: LOAD_MEM.bytes)
        and:
        def trace = Mock(TraceRecord)

        when:
        stats.markCompleted(task,trace)

        then:
        task.aborted >> true
        task.getHashLog() >> HASH
        and:
        stats.runningCount == RUNNING -1
        stats.loadCpus == LOAD_CPUS - CPUS
        stats.loadMemory == (LOAD_MEM - MEM).bytes
        stats.abortedCount == ABORTED +1
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
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
            getConfig() >> Mock(TaskConfig) { getCpus() >> CPUS; getMemory() >> MEM } }
        and:
        def rec = new ProgressRecord(0, 'foo')
        rec.running = RUNNING
        rec.failed = FAILED
        rec.retries = RETRIES
        rec.ignored = IGNORED
        rec.aborted = ABORTED
        rec.succeeded = SUCCEEDED
        rec.loadCpus = LOAD_CPUS
        rec.loadMemory = LOAD_MEM.bytes
        and:
        def stats = new WorkflowStats(
                records: [0:rec],
                runningCount: RUNNING,
                failedCount: FAILED,
                retriesCount: RETRIES,
                ignoredCount: IGNORED,
                succeededCount: SUCCEEDED,
                abortedCount: ABORTED,
                loadCpus: LOAD_CPUS,
                loadMemory: LOAD_MEM.bytes)
        and:
        def trace = Mock(TraceRecord)

        when:
        stats.markCompleted(task,trace)

        then:
        task.aborted >> false
        task.failed >> false
        task.getHashLog() >> HASH
        and:
        stats.runningCount == RUNNING -1
        stats.succeededCount == SUCCEEDED +1
        stats.loadCpus == LOAD_CPUS - CPUS
        stats.loadMemory == (LOAD_MEM - MEM).bytes

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
        def rec = new ProgressRecord(0, 'foo')
        and:
        assert !rec.terminated
        and:
        def stats = new WorkflowStats( records: [0:rec] )
        when:
        stats.markTerminated(Mock(TaskProcessor) { getId() >> 0 } )
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
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
        }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.cached = CACHED
        rec.stored = STORED
        and:
        def trace = Mock(TraceRecord)
        def stats = new WorkflowStats(
                records: [0:rec],
                cachedCount: CACHED )

        when:
        stats.markCached(task, trace)
        then:
        1 * trace.get('realtime')  >> 5000
        1 * trace.get('cpus') >> 1
        and:
        stats.cachedCount == CACHED +1
        stats.cachedDuration == 5.sec
        and:
        rec.hash == 'XYZ'
        rec.cached == CACHED +1
        rec.stored == STORED
    }

    def 'should mark store' () {
        given:
        def CACHED = 10
        def STORED = 20
        and:
        def task = Mock(TaskRun) {
            getHashLog() >> 'XYZ'
            getProcessor() >> Mock(TaskProcessor) { getId() >> 0 }
        }
        and:
        def rec = new ProgressRecord(10, 'foo')
        rec.cached = CACHED
        rec.stored = STORED
        and:
        def trace = Mock(TraceRecord)
        def stats = new WorkflowStats(
                records: [0:rec],
                cachedCount: CACHED )

        when:
        stats.markCached(task, null)
        then:
        stats.cachedCount == CACHED
        and:
        rec.cached == CACHED
        rec.stored == STORED +1
        and:
        rec.hash == 'skipped'
    }


    def 'should check records len' () {
        given:
        def p1 = Mock(TaskProcessor) { getId() >> 1; getName() >> 'foo' }
        def p2 = Mock(TaskProcessor) { getId() >> 5; getName() >> 'bar' }

        when:
        def stats = new WorkflowStats()
        then:
        stats.progressLength == 0

        when:
        stats.markCreated(p1)
        then:
        stats.getProgressLength() == 1

        when:
        stats.markPending(p2)
        then:
        stats.getProgressLength() == 2

    }
}
