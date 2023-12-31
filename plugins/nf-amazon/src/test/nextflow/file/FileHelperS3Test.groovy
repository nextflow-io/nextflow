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

import nextflow.SysEnv
import spock.lang.Unroll
import test.AppSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperS3Test extends AppSpec {

    private Path asPath(String str) {
        if( str==null )
            return null
        if( str.startsWith('s3://'))
            return FileSystemPathFactory.parse(str)
        else
            return Path.of(str)
    }

    @Unroll
    def 'should convert to canonical path with base' () {
        given:
        SysEnv.push(NXF_FILE_ROOT: 's3://host.com/work')

        expect:
        FileHelper.toCanonicalPath(VALUE) == asPath(EXPECTED)

        cleanup:
        SysEnv.pop()

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | 's3://host.com/work/file.txt'
        Path.of('file.txt')         | 's3://host.com/work/file.txt'
        and:
        './file.txt'                | 's3://host.com/work/file.txt'
        '.'                         | 's3://host.com/work'
        './'                        | 's3://host.com/work'
        '../file.txt'               | 's3://host.com/file.txt'
        and:
        '/file.txt'                 | '/file.txt'
        Path.of('/file.txt')        | '/file.txt'

    }

}
