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

package nextflow.executor


import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptOutputFilesTest extends Specification {

    @Unroll
    def 'should escape output file names' () {

        expect:
        ScriptOutputFiles.shellSpecialChars(NAME, ScriptOutputFiles.glob(GLOB)) == EXPECTED

        where:
        NAME        | GLOB      | EXPECTED
        'foo.txt'   | false     | 'foo.txt'
        'foo.txt'   | true      | 'foo.txt'
        and:
        'fo .txt'   | false     | 'fo\\ .txt'
        'fo .txt'   | true      | 'fo\\ .txt'
        and:
        'fo*.txt'   | false     | 'fo\\*.txt'
        'fo*.txt'   | true      | 'fo*.txt'
        'fo?.txt'   | false     | 'fo\\?.txt'
        'fo?.txt'   | true      | 'fo?.txt'
    }

    def 'should return shell escaped names' () {
        given:
        def files = new ScriptOutputFiles()
        and:
        files.putName('foo.txt', true)
        files.putName('foo-*.txt', true)    // the * should be interpreted as a glob => not escaped
        files.putName('foo-?.txt', false)   // the ? should not be interpreted as a glob => escaped
        files.putName('f o o.txt', false)   // glob does not matter => blank should be escaped
        
        expect:
        files.toShellEscapedNames() == 'foo.txt foo-*.txt foo-\\?.txt f\\ o\\ o.txt'

    }

}
