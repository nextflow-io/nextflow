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

package nextflow.script

import java.nio.file.Files
import java.nio.file.Path

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.RecordMap
import nextflow.util.VersionNumber
import spock.lang.Specification

/**
 * Tests for {@link ParamsHelper}.
 *
 * This is the shared mapping of a param value to a declared type, used by
 * the `params` block, typed workflows, and typed processes. The exhaustive
 * per-type coverage lives here -- each caller only needs to prove that it
 * routes through this class.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ParamsHelperTest extends Specification {

    private static Param param(java.lang.reflect.Type type) {
        new Param('x', type, false, null)
    }

    def 'should resolve a value against the declared type'() {
        expect:
        final result = FROM_CLI
            ? ParamsHelper.resolveFromCli(param(TYPE), VALUE)
            : ParamsHelper.resolveFromCode(param(TYPE), VALUE)
        result == EXPECTED
        result?.getClass() == EXPECTED?.getClass()

        where:
        FROM_CLI | TYPE          | VALUE     | EXPECTED
        // a command-line value is always a string, so it is parsed
        true     | String        | 'hello'   | 'hello'
        true     | Boolean       | 'true'    | true
        true     | Boolean       | 'false'   | false
        true     | Boolean       | 'TRUE'    | true
        true     | Integer       | '42'      | 42
        true     | Integer       | '0'       | 0
        true     | Integer       | ' 42 '    | 42
        true     | Float         | '5.5'     | 5.5f
        true     | Float         | '5'       | 5.0f
        true     | Duration      | '2h'      | Duration.of('2h')
        true     | MemoryUnit    | '4.GB'    | MemoryUnit.of('4.GB')
        true     | VersionNumber | '1.2.3'   | new VersionNumber('1.2.3')
        // a config value is already structured, but a number is still
        // normalized to the declared type
        false    | String        | 'hello'   | 'hello'
        false    | Integer       | 40        | 40
        false    | Integer       | 40G       | 40
        false    | Integer       | 40.0G     | 40
        false    | Float         | 96.4G     | 96.4f
        false    | Float         | 96.4d     | 96.4f
        false    | Float         | 40        | 40.0f
        false    | Duration      | '2h'      | Duration.of('2h')
        false    | MemoryUnit    | '4.GB'    | MemoryUnit.of('4.GB')
        false    | VersionNumber | '1.2.3'   | new VersionNumber('1.2.3')
        // a null value is always passed through
        true     | Integer       | null      | null
        false    | Integer       | null      | null
        // an untyped param is passed through unchanged
        true     | null          | '42'      | '42'
    }

    def 'should widen an integral value that is too large for an Integer'() {
        expect:
        final result = ParamsHelper.resolveFromCli(param(Integer), VALUE)
        result == EXPECTED
        result.getClass() == EXPECTED.getClass()

        where:
        VALUE                          | EXPECTED
        '2147483647'                   | 2147483647                       // Integer.MAX_VALUE
        '99999999999'                  | 99999999999L
        '-99999999999'                 | -99999999999L
        '999999999999999999999999999'  | 999999999999999999999999999G
    }

    def 'should leave a value that the declared type cannot represent unconverted'() {
        when:
        // the value is reported by the caller as not assignable to the
        // declared type, rather than being silently coerced
        final result = ParamsHelper.resolveFromCli(param(TYPE), VALUE)

        then:
        result == VALUE
        !ParamsHelper.isAssignableFrom(TYPE, result.getClass())

        where:
        TYPE    | VALUE
        Integer | 'abc'
        Integer | '3.7'      // a fractional value is not truncated
        Float   | 'abc'
    }

    def 'should not convert a boolean string that came from the config'() {
        expect:
        // unlike a command-line value, a config value is expected to already
        // be a boolean -- a string is left alone and reported as a mismatch
        ParamsHelper.resolveFromCode(param(Boolean), 'false') == 'false'
        ParamsHelper.resolveFromCode(param(Boolean), false) == false
    }

    def 'should resolve a Path and require it to exist'() {
        given:
        def file = Files.createTempFile('data', '.txt')

        when:
        def result = ParamsHelper.resolveFromCli(param(Path), file.toString())

        then:
        result == file

        when:
        ParamsHelper.resolveFromCli(param(Path), '/some/missing/file.txt')

        then:
        def e = thrown(Exception)
        e.message.contains('/some/missing/file.txt')

        cleanup:
        file?.delete()
    }

    def 'should determine whether a value type is assignable to a declared type'() {
        expect:
        ParamsHelper.isAssignableFrom(TARGET, SOURCE) == EXPECTED

        where:
        TARGET  | SOURCE     | EXPECTED
        Integer | Integer    | true
        Integer | Long       | true
        Integer | BigInteger | true
        Integer | Float      | false
        Integer | String     | false
        Float   | Float      | true
        Float   | Integer    | true
        Float   | BigDecimal | true
        Float   | String     | false
        String  | String     | true
        String  | Integer    | false
        SampleRec | RecordMap  | true
        SampleRec | Map        | false
    }

    static class SampleRec implements nextflow.script.types.Record {
        String id
    }

}
