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
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperS3Test extends Specification {


    @Unroll
    def 'should convert to canonical path with base' () {
        given:
        SysEnv.push(NXF_FILE_ROOT: 's3://host.com/work')

        expect:
        FileHelper.toCanonicalPath(VALUE) == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | FileSystemPathFactory.parse('s3://host.com/work/file.txt')
        Path.of('file.txt')         | FileSystemPathFactory.parse('s3://host.com/work/file.txt')
        and:
        './file.txt'                | FileSystemPathFactory.parse('s3://host.com/work/file.txt')
        '.'                         | FileSystemPathFactory.parse('s3://host.com/work')
        './'                        | FileSystemPathFactory.parse('s3://host.com/work')
        '../file.txt'               | FileSystemPathFactory.parse('s3://host.com/file.txt')
        and:
        '/file.txt'                 | Path.of('/file.txt')
        Path.of('/file.txt')        | Path.of('/file.txt')

    }

    def 'should convert to a canonical path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }

        expect:
        FileHelper.toCanonicalPath(VALUE).toUri() == EXPECTED

        where:
        VALUE                       | EXPECTED
        's3://foo/some/file.txt'    | new URI('s3:///foo/some/file.txt')
        's3://foo/some///file.txt'  | new URI('s3:///foo/some/file.txt')
    }

    @Unroll
    def 'should remove consecutive slashes in the path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }

        expect:
        FileHelper.asPath(STR).toUri() == EXPECTED
        where:
        STR                         | EXPECTED
        's3://foo//this/that'       | new URI('s3:///foo/this/that')
        's3://foo//this///that'     | new URI('s3:///foo/this/that')
    }
}
