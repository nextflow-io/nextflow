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
import spock.lang.Unroll
import test.AppSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class FileHelperAzTest extends AppSpec {

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

    private Path asPath(String str) {
        if( str==null )
            return null
        if( str.startsWith('az://'))
            return FileSystemPathFactory.parse(str)
        if( str.startsWith('/'))
            return Path.of(str)
        throw new IllegalStateException("Unexpected file path: $str")
    }

    @Unroll
    def 'should convert to canonical path with base' () {
        given:
        SysEnv.push(NXF_FILE_ROOT: 'az://host.com/work')

        expect:
        FileHelper.toCanonicalPath(VALUE) == asPath(EXPECTED)

        cleanup:
        SysEnv.pop()

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | 'az://host.com/work/file.txt'
        Path.of('file.txt')         | 'az://host.com/work/file.txt'
        and:
        './file.txt'                | 'az://host.com/work/file.txt'
        '.'                         | 'az://host.com/work'
        './'                        | 'az://host.com/work'
        '../file.txt'               | 'az://host.com/file.txt'
        and:
        '/file.txt'                 | '/file.txt'
        Path.of('/file.txt')        | '/file.txt'

    }

}
