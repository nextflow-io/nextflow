/*
 * Copyright 2013-2026, Seqera Labs
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

import nextflow.Session
import nextflow.SysEnv
import nextflow.trace.event.TaskEvent
import spock.lang.Specification
import spock.lang.Unroll

class TerminalActivityTrackerTest extends Specification {

    static final char ESC = 0x1B as char

    def 'should emit only when the activity state changes'() {
        given:
        def sink = new StringBuilder()
        def tracker = TerminalActivityTracker.open(new Osc7501Emitter(sink))

        when:
        tracker.onStats(0)
        tracker.onStats(0)
        tracker.onStats(4)
        tracker.onStats(2)
        tracker.onStats(1)
        tracker.onStats(0)
        tracker.onStats(3)

        then:
        sink.toString() == sequence('idle', 'Waiting for tasks') +
                sequence('working', 'Running tasks') +
                sequence('idle', 'Waiting for tasks') +
                sequence('working', 'Running tasks')
        sink.toString().count('state=working') == 2
        sink.toString().count('state=idle') == 2
    }

    def 'should report completion after the queue drains without an idle flicker'() {
        given:
        def sink = new StringBuilder()
        def tracker = TerminalActivityTracker.open(new Osc7501Emitter(sink))
        tracker.onStats(2)

        when:
        tracker.beginShutdown()
        tracker.onStats(0)
        tracker.onRunComplete(true, false)

        then:
        sink.toString() == sequence('working', 'Running tasks') +
                sequence('done', 'Workflow completed')
    }

    def 'should report a stopped failure once'() {
        given:
        def sink = new StringBuilder()
        def tracker = TerminalActivityTracker.open(new Osc7501Emitter(sink))
        tracker.onStats(1)

        when:
        tracker.onFlowError(false)
        tracker.onRunComplete(false, false)
        tracker.onStats(1)

        then:
        sink.toString() == sequence('working', 'Running tasks') +
                sequence('error', 'Workflow failed')
    }

    def 'should report a cancel as idle'() {
        given:
        def sink = new StringBuilder()
        def tracker = TerminalActivityTracker.open(new Osc7501Emitter(sink))
        tracker.onStats(1)

        when:
        tracker.onFlowError(true)

        then:
        sink.toString() == sequence('working', 'Running tasks') +
                sequence('idle', 'Run cancelled')
    }

    def 'should report a failed completion when no task error arrived first'() {
        given:
        def sink = new StringBuilder()
        def tracker = TerminalActivityTracker.open(new Osc7501Emitter(sink))

        when:
        tracker.onStats(1)
        tracker.beginShutdown()
        tracker.onRunComplete(false, false)

        then:
        sink.toString().endsWith(sequence('error', 'Workflow failed'))
    }

    @Unroll
    def 'should enable reporting only for an ANSI tty that has not opted out'() {
        given:
        SysEnv.push(ENV as Map<String,String>)
        def session = Stub(Session) {
            isAnsiLog() >> ANSI
            getConfig() >> CONFIG
        }

        expect:
        TerminalActivityTracker.enabled(session, TTY) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        ENV                      | ANSI  | TTY   | CONFIG            | EXPECTED
        [:]                      | true  | true  | [:]               | true
        [:]                      | true  | true  | null              | true
        [:]                      | false | true  | [:]               | false
        [:]                      | true  | false | [:]               | false
        [NXF_OSC7501: 'false']   | true  | true  | [:]               | false
        [NXF_OSC7501: '0']       | true  | true  | [:]               | false
        [NXF_OSC7501: 'true']    | true  | true  | [:]               | true
        [NXF_OSC7501: '1']       | true  | true  | [:]               | true
        [NXF_OSC7501: '']        | true  | true  | [:]               | true
        [:]                      | true  | true  | [osc7501: false]  | false
        [:]                      | true  | true  | [osc7501: 'false']| false
        [:]                      | true  | true  | [osc7501: true]   | true
        [NXF_OSC7501: 'false']   | true  | true  | [osc7501: true]   | false
        [NXF_OSC7501: 'true']    | true  | true  | [osc7501: false]  | false
    }

    def 'should forward a workflow failure and a clean completion as terminal bytes'() {
        given:
        def failed = new StringBuilder()
        def failedObserver = observer(failed)
        failedObserver.@session = Stub(Session) {
            isCancelled() >> false
            getError() >> new RuntimeException('task failed')
            isSuccess() >> false
        }
        failedObserver.publishActivity(stats(2))

        when:
        failedObserver.onFlowError(Stub(TaskEvent))
        failedObserver.finishActivity()

        then:
        failed.toString() == sequence('working', 'Running tasks') +
                sequence('error', 'Workflow failed')
        and:
        !AnsiLogObserver.userInterrupt(failedObserver.@session as Session)

        when:
        def done = new StringBuilder()
        def doneObserver = observer(done)
        doneObserver.@session = Stub(Session) {
            isCancelled() >> false
            getError() >> null
            isSuccess() >> true
        }
        doneObserver.publishActivity(stats(1))
        doneObserver.@activity.beginShutdown()
        doneObserver.publishActivity(stats(0))
        doneObserver.finishActivity()

        then:
        done.toString() == sequence('working', 'Running tasks') +
                sequence('done', 'Workflow completed')
    }

    def 'should treat a cancelled session as a user interrupt'() {
        given:
        def session = Stub(Session) {
            isCancelled() >> true
            getError() >> new RuntimeException('ignored')
        }
        def sink = new StringBuilder()
        def observer = observer(sink)
        observer.@session = session
        observer.publishActivity(stats(1))

        when:
        observer.onFlowError(Stub(TaskEvent))

        then:
        AnsiLogObserver.userInterrupt(session)
        sink.toString() == sequence('working', 'Running tasks') +
                sequence('idle', 'Run cancelled')
    }

    private static AnsiLogObserver observer(StringBuilder sink) {
        def observer = new AnsiLogObserver()
        observer.@activity = TerminalActivityTracker.open(new Osc7501Emitter(sink))
        return observer
    }

    private static WorkflowStats stats(int running) {
        def stats = new WorkflowStats()
        stats.@runningCount = running
        return stats
    }

    private static String sequence(String state, String message) {
        return Osc7501Emitter.format(state, message)
    }

}
