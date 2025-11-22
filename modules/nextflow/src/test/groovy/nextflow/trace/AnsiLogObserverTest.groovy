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

}
