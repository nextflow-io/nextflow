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

package io.seqera.wave.plugin.cli

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveDebugCmdTest extends Specification {

    @Unroll
    def 'should check remote path' () {
        expect:
        WaveDebugCmd.isRemotePath(PATH) == EXPECTED
        where:
        PATH                | EXPECTED
        null                | false
        'foo'               | false
        '/some/file'        | false
        and:
        's3://foo/bar'      | true
        'gs://foo/bar'      | true
        and:
        'file:/foo/bar'     | false
        'file://foo/bar'    | false
        'file:///foo/bar'   | false
    }

}
