/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.file

import java.nio.file.Path

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class FileHelperAzTest extends Specification {

    def setupSpec() {
        def CONFIG = [azure: [
                storage: [
                        accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                        accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
                ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
    }

    def cleanupSpec() {
        Global.session = null
    }


    @Unroll
    def 'should convert to canonical path with base' () {
        given:
        SysEnv.push(NXF_FILE_ROOT: 'az://host.com/work')

        expect:
        FileHelper.toCanonicalPath(VALUE) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | FileHelper.asPath('az://host.com/work/file.txt')
        Path.of('file.txt')         | FileHelper.asPath('az://host.com/work/file.txt')
        and:
        './file.txt'                | FileHelper.asPath('az://host.com/work/file.txt')
        '.'                         | FileHelper.asPath('az://host.com/work')
        './'                        | FileHelper.asPath('az://host.com/work')
        '../file.txt'               | FileHelper.asPath('az://host.com/file.txt')
        and:
        '/file.txt'                 | Path.of('/file.txt')
        Path.of('/file.txt')        | Path.of('/file.txt')

    }

    def 'should convert to a canonical path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [google:[project:'foo', region:'x']] }

        expect:
        FileHelper.toCanonicalPath(VALUE).toUri() == EXPECTED

        where:
        VALUE                       | EXPECTED
        'az://foo/some/file.txt'    | new URI('az://foo/some/file.txt')
        'az://foo/some///file.txt'  | new URI('az://foo/some/file.txt')
    }

    @Unroll
    def 'should remove consecutive slashes in the path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }

        expect:
        FileHelper.asPath(STR).toUri() == EXPECTED
        where:
        STR                         | EXPECTED
        'az://foo//this/that'       | new URI('az://foo/this/that')
        'az://foo//this///that'     | new URI('az://foo/this/that')
    }
}
