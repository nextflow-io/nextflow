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

package nextflow.trace

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import nextflow.trace.event.TaskEvent
import spock.lang.Specification

/**
 * Tests for AgentLogObserver
 *
 * @author Edmund Miller <edmund.miller@utdallas.edu>
 */
class AgentLogObserverTest extends Specification {

    def 'should output pipeline info on flow begin'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)
        def session = Stub(Session) {
            getManifest() >> [name: 'nf-core/rnaseq', version: '3.14.0']
            getScriptName() >> 'main.nf'
            getProfile() >> 'test,docker'
            getWorkDir() >> Paths.get('/tmp/work')
        }

        when:
        observer.onFlowCreate(session)
        observer.onFlowBegin()

        then:
        output[0] == '[PIPELINE] nf-core/rnaseq 3.14.0 | profile=test,docker'
        output[1].startsWith('[WORKDIR]')
        output[1].contains('/tmp/work')
    }

    def 'should deduplicate warnings'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)

        when:
        observer.appendWarning('Task runtime metrics are not reported when using macOS')
        observer.appendWarning('Task runtime metrics are not reported when using macOS')
        observer.appendWarning('Different warning message')

        then:
        output.size() == 2
        output[0] == '[WARN] Task runtime metrics are not reported when using macOS'
        output[1] == '[WARN] Different warning message'
    }

    def 'should output error with context'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)
        def task = Stub(TaskRun) {
            getName() >> 'FASTQC (sample_1)'
            getExitStatus() >> 127
            getScript() >> 'fastqc --version'
            getStderr() >> 'bash: fastqc: command not found'
            getStdout() >> null
            getWorkDir() >> Paths.get('/tmp/work/ab/123456')
        }

        when:
        observer.printTaskError(task)

        then:
        output[0] == '[ERROR] FASTQC (sample_1)'
        output[1] == 'exit: 127'
        output[2] == 'cmd: fastqc --version'
        output[3] == 'stderr: bash: fastqc: command not found'
        output[4].startsWith('workdir:')
        output[4].contains('/tmp/work/ab/123456')
    }

    def 'should output summary on flow complete'() {
        given:
        def output = []
        def stats = new WorkflowStats()
        stats.succeededCount = 10
        stats.failedCount = 0
        stats.cachedCount = 5
        def statsObserver = Stub(WorkflowStatsObserver) {
            getStats() >> stats
        }
        def observer = new TestAgentLogObserver(output)
        observer.setStatsObserver(statsObserver)

        when:
        observer.onFlowComplete()

        then:
        output[0] == '\n[SUCCESS] completed=10 failed=0 cached=5'
    }

    def 'should output failed summary when tasks fail'() {
        given:
        def output = []
        def stats = new WorkflowStats()
        stats.succeededCount = 5
        stats.failedCount = 2
        stats.cachedCount = 0
        def statsObserver = Stub(WorkflowStatsObserver) {
            getStats() >> stats
        }
        def observer = new TestAgentLogObserver(output)
        observer.setStatsObserver(statsObserver)

        when:
        observer.onFlowComplete()

        then:
        output[0] == '\n[FAILED] completed=7 failed=2 cached=0'
    }

    def 'should handle task complete for failed task'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)
        def task = Stub(TaskRun) {
            getName() >> 'TEST_PROC'
            getExitStatus() >> 1
            getScript() >> 'exit 1'
            getStderr() >> 'error'
            getStdout() >> null
            getWorkDir() >> Paths.get('/work/xx/yy')
            isFailed() >> true
        }
        def handler = Stub(TaskHandler) {
            getTask() >> task
        }
        def event = Stub(TaskEvent) {
            getHandler() >> handler
        }

        when:
        observer.onTaskComplete(event)

        then:
        output.size() > 0
        output[0] == '[ERROR] TEST_PROC'
    }

    def 'should not output for successful task'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)
        def task = Stub(TaskRun) {
            isFailed() >> false
        }
        def handler = Stub(TaskHandler) {
            getTask() >> task
        }
        def event = Stub(TaskEvent) {
            getHandler() >> handler
        }

        when:
        observer.onTaskComplete(event)

        then:
        output.size() == 0
    }

    def 'should truncate long command'() {
        given:
        def output = []
        def observer = new TestAgentLogObserver(output)
        def longCommand = 'x' * 300
        def task = Stub(TaskRun) {
            getName() >> 'PROC'
            getExitStatus() >> 1
            getScript() >> longCommand
            getStderr() >> null
            getStdout() >> null
            getWorkDir() >> Paths.get('/work')
        }

        when:
        observer.printTaskError(task)

        then:
        def cmdLine = output.find { it.startsWith('cmd:') }
        cmdLine != null
        cmdLine.length() < 250
        cmdLine.endsWith('...')
    }

    /**
     * Test subclass that captures output
     */
    static class TestAgentLogObserver extends AgentLogObserver {
        private List<String> output

        TestAgentLogObserver(List<String> output) {
            this.output = output
        }

        @Override
        protected void printLine(String line) {
            output << line
        }
    }
}
