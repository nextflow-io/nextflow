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
        CmdLineHelper.parseGnuArgs(SOURCE).toString()  == EXPECTED

        where:
        SOURCE              | EXPECTED
        'a b c'             | '[]'
        'a -b c'            | '[option{b: [c]}]'
        'a -b -1'           | '[option{b: [-1]}]'
        '-a b -1'           | '[option{a: [b, -1]}]'
        "-a 'b -1'"         | '[option{a: [b -1]}]'
        "-a='b -1'"         | '[option{a: [b -1]}]'
        '-a "b c"'          | '[option{a: [b c]}]'
        '-a="b c"'          | '[option{a: [b c]}]'
        and:
        '--foo 1'               | '[option{foo: [1]}]'
        '--foo 1 --bar 2'       | '[option{foo: [1]}, option{bar: [2]}]'
        '--foo 1 2 3 --bar 4'   | '[option{foo: [1, 2, 3]}, option{bar: [4]}]'
        '--foo 1 2 3 --bar'     | '[option{foo: [1, 2, 3]}, option{bar: [true]}]'
        '--foo --bar'           | '[option{foo: [true]}, option{bar: [true]}]'
        and:
        // single non-gnu is not capture
        '--foo 1 -bar 2'        | '[option{foo: [1, -bar, 2]}]'
        and:
        '--foo-name 1 --bar-opt 2' | '[option{foo-name: [1]}, option{bar-opt: [2]}]'
    }

}
