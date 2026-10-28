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

package nextflow.cli

import spock.lang.Specification

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConsoleInputTest extends Specification {

    def 'should read all the lines from the standard input'() {
        given:
        def stdin = System.in
        System.setIn(new ByteArrayInputStream('alpha\nbeta\ngamma\n'.bytes))
        and:
        def input = new ConsoleInput()

        expect:
        input.readLine() == 'alpha'
        input.readLine() == 'beta'
        input.readLine() == 'gamma'
        input.readLine() == null

        cleanup:
        System.setIn(stdin)
    }

    def 'should return null when the standard input is empty'() {
        given:
        def stdin = System.in
        System.setIn(new ByteArrayInputStream(new byte[0]))

        expect:
        new ConsoleInput().readLine() == null

        cleanup:
        System.setIn(stdin)
    }

}
