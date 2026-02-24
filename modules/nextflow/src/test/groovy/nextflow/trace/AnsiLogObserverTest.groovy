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

import nextflow.Session
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

    def 'should create hyperlink' () {
        expect:
        // Local paths (starting with /) get file:// prefix added automatically
        AnsiLogObserver.hyperlink('hash', '/path/to/work') == '\033]8;;file:///path/to/work\007hash\033]8;;\007'
        // URLs with schemes are used as-is
        AnsiLogObserver.hyperlink('hash', 's3://bucket/path') == '\033]8;;s3://bucket/path\007hash\033]8;;\007'
        AnsiLogObserver.hyperlink('hash', 'gs://bucket/path') == '\033]8;;gs://bucket/path\007hash\033]8;;\007'
        AnsiLogObserver.hyperlink('hash', 'az://container/path') == '\033]8;;az://container/path\007hash\033]8;;\007'
        AnsiLogObserver.hyperlink('text', null) == 'text'
        AnsiLogObserver.hyperlink('text', '') == 'text'
    }

    def 'should render hash as hyperlink when workDir is set' () {
        given:
        def session = Mock(Session) { getConfig() >> [cleanup: false] }
        def observer = new AnsiLogObserver()
        observer.@session = session
        observer.@labelWidth = 3
        observer.@cols = 190
        and:
        def stats = new ProgressRecord(1, 'foo')
        stats.submitted = 1
        stats.hash = '4e/486876'
        stats.workDir = WORKDIR

        when:
        def result = observer.line(stats).toString()

        then:
        result.contains('\033]8;;' + EXPECTED_HREF + '\007')
        result.contains('\033]8;;\007')

        where:
        WORKDIR                          | EXPECTED_HREF
        '/work/4e/486876abc'             | 'file:///work/4e/486876abc'
        's3://bucket/work/4e/486876abc'  | 's3://bucket/work/4e/486876abc'
    }

    def 'should strip ansi escape codes' () {
        expect:
        AnsiLogObserver.stripAnsi('hello') == 'hello'
        AnsiLogObserver.stripAnsi('\033[32mgreen\033[0m') == 'green'
        AnsiLogObserver.stripAnsi('\033[1;31mbold red\033[0m') == 'bold red'
        AnsiLogObserver.stripAnsi('\033]8;;http://example.com\007link\033]8;;\007') == 'link'
        AnsiLogObserver.stripAnsi('\033[2m[\033[0m\033[34mab/123456\033[0m\033[2m] \033[0mfoo') == '[ab/123456] foo'
    }

    @Unroll
    def 'should count visual lines with wrapping' () {
        given:
        def observer = new AnsiLogObserver()
        observer.@cols = COLS

        expect:
        observer.countVisualLines(INPUT) == EXPECTED

        where:
        COLS | INPUT                        | EXPECTED
        80   | 'short line\n'               | 1
        80   | 'line1\nline2\n'             | 2
        80   | 'a' * 80 + '\n'              | 1       // exactly fits, no wrap
        80   | 'a' * 81 + '\n'              | 2       // wraps to 2 lines
        80   | 'a' * 160 + '\n'             | 2       // exactly 2 lines
        80   | 'a' * 161 + '\n'             | 3       // wraps to 3 lines
        40   | 'a' * 100 + '\n'             | 3       // 100 chars in 40-col terminal
        80   | 'short\n' + 'a' * 200 + '\n'| 4       // 1 + 3 lines
    }

    def 'should not render hyperlink when cleanup is enabled' () {
        given:
        def session = Mock(Session) { getConfig() >> [cleanup: true] }
        def observer = new AnsiLogObserver()
        observer.@session = session
        observer.@labelWidth = 3
        observer.@cols = 190
        and:
        def stats = new ProgressRecord(1, 'foo')
        stats.submitted = 1
        stats.hash = '4e/486876'
        stats.workDir = '/work/4e/486876abc'

        when:
        def result = observer.line(stats).toString()

        then:
        // Should NOT contain hyperlink start sequence
        !result.contains('\033]8;;')
    }

}
