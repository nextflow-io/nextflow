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

    def 'should not implement a path matcher' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        def path = FileHelper.asPath('s3://my-bucket/work/ab/cdef/.fusion')

        when:
        path.getFileSystem().getPathMatcher('glob:*')
        then:
        // this is the reason why the `getDefaultPathMatcher` fallback is needed
        thrown(UnsupportedOperationException)
    }

    def 'should return a bucket-less file name' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        def path = FileHelper.asPath('s3://my-bucket/work/ab/cdef/.fusion')

        expect:
        // the name must not be turned into a bucket, otherwise it could not be matched
        path.getFileName().toString() == '.fusion'
        and:
        !path.getFileName().isAbsolute()
        and:
        path.getFileSystem().getPath('.fusion').toString() == '.fusion'
    }

    @Unroll
    def 'should match the hidden component #COMPONENT of a s3 path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [:] }
        def path = FileHelper.asPath(PATH)
        and:
        // the S3 file system does not implement `getPathMatcher`, hence this is the
        // `getDefaultPathMatcher` fallback, matching against the path `toString()`
        def matchers = [ FileHelper.getPathMatcherFor("glob:$COMPONENT", path.getFileSystem()) ]

        expect:
        FileHelper.matchesAnyName(matchers, path) == EXPECTED

        where:
        COMPONENT   | PATH                                          | EXPECTED
        '.fusion'   | 's3://my-bucket/work/ab/cdef/.fusion'         | true
        '.fusion'   | 's3://my-bucket/.fusion'                      | true
        '.fusion'   | 's3://my-bucket/work/ab/cdef/.fusion-other'   | false
        '.fusion'   | 's3://my-bucket/work/ab/cdef'                 | false
        '.fusion'   | 's3://.fusion/work/ab/cdef'                   | false
        '.*'        | 's3://my-bucket/work/ab/cdef/.fusion'         | true
        '.*'        | 's3://my-bucket/work/ab/.git'                 | true
        '\\.bar'    | 's3://my-bucket/work/ab/.bar'                 | true
        '.b*'       | 's3://my-bucket/work/ab/.bar'                 | true
        '.b*'       | 's3://my-bucket/work/ab/.zzz'                 | false
    }
}
