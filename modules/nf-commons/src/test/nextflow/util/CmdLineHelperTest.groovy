/*
 * Copyright 2020-2021, Seqera Labs
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
 *
 */

package nextflow.util

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdLineHelperTest extends Specification{

    @Unroll
    def 'should parse command line' () {
        expect:
        CmdLineHelper.splitter(SOURCE)  == EXPECTED

        where:
        SOURCE              | EXPECTED
        "foo bar baz"       | ['foo','bar','baz']
        "foo 'this that'"   | ['foo', 'this that']
    }

    @Unroll
    def 'should parse args' () {
        expect:
        CmdLineHelper.parseGnuArgs(SOURCE)  == EXPECTED

        where:
        SOURCE              | EXPECTED
        'a b c'             | [:]
        'a -b c'            | [b:'c']
        'a -b -1'           | [b:'-1']
        '-a b -1'           | [a:'b -1']
        "-a 'b -1'"         | [a:'b -1']
        "-a='b -1'"         | [a:'b -1']
        '-a "b c"'          | [a:'b c']
        '-a="b c"'          | [a:'b c']
        and:
        '--foo 1'               | [foo: '1']
        '--foo 1 --bar 2'       | [foo: '1', bar: '2']
        '--foo 1 2 3 --bar 4'   | [foo: '1 2 3', bar: '4']
        '--foo 1 2 3 --bar'     | [foo: '1 2 3', bar: 'true']
        '--foo --bar'           | [foo: 'true', bar: 'true']
        and:
        // single non-gnu is not capture
        '--foo 1 -bar 2'        | [foo: '1 -bar 2']
    }

}
