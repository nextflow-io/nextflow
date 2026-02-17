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

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AnsiLogObserverTest extends Specification {

    @Unroll
    def 'should render a process line' () {

        given:
        def observer = new AnsiLogObserver()
        def stats = new ProgressRecord(1, 'foo')
        stats.submitted = SUBMIT
        stats.succeeded = SUCCEEDED
        stats.cached    = CACHE
        stats.stored    = STORE
        stats.terminated = DONE
        stats.errored = ERR
        stats.hash = HASH

        when:
        observer.@labelWidth = stats.name.size()
        observer.@cols = WIDTH

        then:
        observer.line(stats).toString().replaceAll('\u001B\\[[\\d;]*[^\\d;]','') == EXPECTED

        where:
        HASH        | SUBMIT  | SUCCEEDED | CACHE | STORE | DONE  | ERR   | WIDTH | EXPECTED
        null        | 0       | 0         | 0     | 0     | false | false | 190   | '[-        ] process > foo -'
        '4e/486876' | 1       | 0         | 0     | 0     | false | false | 190   | '[4e/486876] process > foo [  0%] 0 of 1'
        '4e/486876' | 0       | 1         | 0     | 0     | false | false | 190   | '[4e/486876] process > foo [100%] 1 of 1'
        '4e/486876' | 1       | 1         | 0     | 0     | false | false | 190   | '[4e/486876] process > foo [ 50%] 1 of 2'
        '4e/486876' | 5       | 5         | 0     | 0     | false | false | 190   | '[4e/486876] process > foo [ 50%] 5 of 10'
        '4e/486876' | 6       | 6         | 0     | 0     | false | false | 180   | '[4e/486876] foo [ 50%] 6 of 12'
        '4e/486876' | 7       | 7         | 0     | 0     | false | false | 70    | '[4e/486876] foo | 7 of 14'
        '4e/486876' | 0       | 0         | 5     | 0     | false | false | 190   | '[4e/486876] process > foo [100%] 5 of 5, cached: 5'
        '4e/486876' | 1       | 1         | 3     | 0     | false | false | 190   | '[4e/486876] process > foo [ 80%] 4 of 5, cached: 3'
        'skipped'   | 0       | 0         | 0     | 5     | false | false | 190   | '[skipped  ] process > foo [100%] 5 of 5, stored: 5'
        'skipped'   | 1       | 1         | 0     | 3     | false | false | 190   | '[skipped  ] process > foo [ 80%] 4 of 5, stored: 3'
        'ab/123456' | 0       | 2         | 0     | 0     | true  | false | 190   | '[ab/123456] process > foo [100%] 2 of 2 ✔'
        'ef/987654' | 0       | 2         | 0     | 0     | true  | true  | 190   | '[ef/987654] process > foo [100%] 2 of 2 ✘'

    }

    def 'should format name' () {
        given:
        def ansi = new AnsiLogObserver()

        expect:
        ansi.fmtWidth(NAME, WIDTH, MAX) == EXPECTED

        where:
        NAME        | WIDTH | MAX  | EXPECTED

        'foo'       | 3     | 80   | 'foo'
        'foo'       | 5     | 80   | 'foo  '

        'long_name' | 9     | 5   | 'long_'
        'long_name' | 9     | 6   | 'lon…me'
        'long_name' | 9     | 2   | 'lo'
        'long_name' | 9     | 3   | 'lon'
        'xx'        | 9     | 1   | 'x'
        'xx'        | 9     | 5   | 'xx   '
        'abcd'      | 9     | 5   | 'abcd '
        '12345678'  | 9     | 5   | '12345'
        '12345678'  | 9     | 6   | '123…78'
    }

    def 'should chop a string' () {
        given:
        def ansi = new AnsiLogObserver()

        expect:
        ansi.fmtChop(NAME, COLS) == EXPECTED

        where:
        NAME        | COLS  | EXPECTED
        'long_name' | 6     | 'lon…me'
        'long_name' | 2     | 'lo'
        'long_name' | 3     | 'lon'
        'xx'        | 1     | 'x'
        'xx'        | 5     | 'xx'
        'abcd'      | 5     | 'abcd'
        '12345678'  | 5     | '12345'
        '12345678'  | 6     | '123…78'

    }

    // -- Terminal notification tests --

    def 'should generate notification sequences' () {
        expect:
        AnsiLogObserver.oscNotify777('Nf', 'done', false) == "\033]777;notify;Nf;done\007"
        AnsiLogObserver.oscNotify777('Nf', 'done', true)  == "\033Ptmux;\033\033]777;notify;Nf;done\007\033\\"
        AnsiLogObserver.oscNotify9('msg', false)           == "\033]9;msg\007"
        AnsiLogObserver.oscNotify9('msg', true)            == "\033Ptmux;\033\033]9;msg\007\033\\"
        AnsiLogObserver.wrapTmuxPassthrough("\033]9;x\007") == "\033Ptmux;\033\033]9;x\007\033\\"
    }

    def 'should generate OSC 99 notification for Kitty' () {
        when:
        def plain = AnsiLogObserver.oscNotify99('T', 'B', false)
        def tmux  = AnsiLogObserver.oscNotify99('T', 'B', true)

        then:
        plain.contains("\033]99;i=nxf:d=0;T\033\\")
        plain.contains("\033]99;i=nxf:p=body;B\033\\")
        tmux.contains("\033Ptmux;")
        tmux.contains("99;i=nxf:d=0;T")
        tmux.contains("99;i=nxf:p=body;B")
    }

    @Unroll
    def 'should detect notify protocol: #DESC' () {
        given:
        SysEnv.push(ENV)

        expect:
        AnsiLogObserver.detectNotifyProtocol() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        DESC                    | ENV                                                      | EXPECTED
        'Kitty'                 | [KITTY_WINDOW_ID: '1']                                   | AnsiLogObserver.NOTIFY_OSC99
        'Kitty over iTerm2'     | [KITTY_WINDOW_ID: '1', TERM_PROGRAM: 'iTerm.app']        | AnsiLogObserver.NOTIFY_OSC99
        'iTerm2 TERM_PROGRAM'   | [TERM_PROGRAM: 'iTerm.app']                              | AnsiLogObserver.NOTIFY_OSC9
        'iTerm2 session ID'     | [ITERM_SESSION_ID: 'w0t0p0:12345']                       | AnsiLogObserver.NOTIFY_OSC9
        'default OSC 777'       | [:]                                                      | AnsiLogObserver.NOTIFY_OSC777
    }

    @Unroll
    def 'should configure from env: #DESC' () {
        given:
        SysEnv.push(ENV)
        def ansi = new AnsiLogObserver()

        expect:
        ansi.@enableOscNotify == NOTIFY
        ansi.@insideTmux == TMUX

        cleanup:
        SysEnv.pop()

        where:
        DESC                | ENV                                              | NOTIFY | TMUX
        'defaults'          | [:]                                              | true   | false
        'notify disabled'   | [NXF_OSC_NOTIFY: 'false']                        | false  | false
        'tmux via TMUX'     | [TMUX: '/tmp/tmux-501/default,12345,0']          | true   | true
        'tmux via TERM'     | [TERM_PROGRAM: 'tmux']                           | true   | true
    }

    @Unroll
    def 'should compute notification body: #DESC' () {
        given:
        SysEnv.push([:])
        def ansi = new AnsiLogObserver()
        ansi.@startTimestamp = START
        ansi.@endTimestamp = END
        ansi.@session = Mock(nextflow.Session) { isSuccess() >> SUCCESS }
        def stats = Mock(WorkflowStats) {
            getSucceedCountFmt() >> SUCCEED
            getCachedCount() >> CACHED
            getCachedCountFmt() >> CACHED_FMT
            getFailedCount() >> FAILED
            getFailedCountFmt() >> FAILED_FMT
        }

        when:
        def result = ansi.computeNotification(stats)

        then:
        result[0] == TITLE
        BODY_CONTAINS.every { result[1].contains(it) }
        BODY_EXCLUDES.every { !result[1].contains(it) }

        cleanup:
        SysEnv.pop()

        where:
        DESC                | SUCCESS | START | END     | SUCCEED | CACHED | CACHED_FMT | FAILED | FAILED_FMT | TITLE                          | BODY_CONTAINS                          | BODY_EXCLUDES
        'success basic'     | true    | 1000  | 61_000  | '10'    | 0      | '0'        | 0      | '0'        | 'Nextflow: completed \u2714'   | ['10 succeeded', '\u2014']             | ['cached', 'failed']
        'failure + counts'  | false   | 1000  | 61_000  | '8'     | 2      | '2'        | 3      | '3'        | 'Nextflow: failed \u2718'      | ['8 succeeded', '2 cached', '3 failed']| []
        'cached only'       | true    | 0     | 5_000   | '0'     | 10     | '10'       | 0      | '0'        | 'Nextflow: completed \u2714'   | ['0 succeeded', '10 cached']           | ['failed']
        'with duration'     | true    | 0     | 150_000 | '5'     | 0      | '0'        | 0      | '0'        | 'Nextflow: completed \u2714'   | ['2m 30s', '\u2014']                   | ['cached', 'failed']
    }

    @Unroll
    def 'should compute notification sequence: #DESC' () {
        given:
        SysEnv.push(ENV)
        def ansi = new AnsiLogObserver()
        ansi.@startTimestamp = 0
        ansi.@endTimestamp = 60_000
        ansi.@session = Mock(nextflow.Session) { isSuccess() >> SUCCESS }
        def stats = Mock(WorkflowStats) {
            getProgressLength() >> 1
            getSucceedCountFmt() >> '5'
            getCachedCount() >> 0
            getFailedCount() >> FAILED
            getFailedCountFmt() >> FAILED.toString()
        }

        when:
        def seq = ansi.computeNotificationSequence(stats)

        then:
        SEQ_CONTAINS.every { seq.contains(it) }

        cleanup:
        SysEnv.pop()

        where:
        DESC                    | ENV                                                      | SUCCESS | FAILED | SEQ_CONTAINS
        'OSC 777 default'       | [:]                                                      | true    | 0      | ['\033]777;notify;', 'completed', '5 succeeded']
        'OSC 777 failure'       | [:]                                                      | false   | 2      | ['\033]777;notify;', 'failed \u2718', '2 failed']
        'OSC 9 iTerm2'          | [TERM_PROGRAM: 'iTerm.app']                              | true    | 0      | ['\033]9;', 'completed \u2714: ']
        'OSC 99 Kitty'          | [KITTY_WINDOW_ID: '1']                                   | true    | 0      | ['99;i=nxf:d=0;', '99;i=nxf:p=body;']
        'OSC 777 + tmux'        | [TMUX: '/tmp/tmux-501/default,1,0']                      | true    | 0      | ['\033Ptmux;', '777;notify;']
        'OSC 9 + tmux'          | [TERM_PROGRAM: 'iTerm.app', TMUX: '/tmp/tmux-501/d,1,0'] | true    | 0      | ['\033Ptmux;', '\033]9;']
        'OSC 99 + tmux'         | [KITTY_WINDOW_ID: '1', TMUX: '/tmp/tmux-501/d,1,0']     | true    | 0      | ['\033Ptmux;', '99;i=nxf:d=0;']
    }

    def 'should return null notification for missing stats' () {
        given:
        SysEnv.push([:])
        def ansi = new AnsiLogObserver()
        def emptyStats = Mock(WorkflowStats) { getProgressLength() >> 0 }

        expect:
        ansi.computeNotificationSequence(null) == null
        ansi.computeNotificationSequence(emptyStats) == null

        cleanup:
        SysEnv.pop()
    }

}
