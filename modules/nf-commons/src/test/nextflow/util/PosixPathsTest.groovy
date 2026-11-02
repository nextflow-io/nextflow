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

package nextflow.util

import java.nio.file.Path

import spock.lang.Specification
import spock.lang.Unroll

class PosixPathsTest extends Specification {

    @Unroll
    def 'should map host path #PATH to #EXPECTED' () {
        expect:
        PosixPaths.toPosix(PATH) == EXPECTED

        where:
        PATH                                | EXPECTED
        null                                | null
        '/home/user/work'                   | '/home/user/work'
        'relative/file.txt'                 | 'relative/file.txt'
        'C:\\Users\\me\\work'               | '/c/Users/me/work'
        'c:\\Users\\me\\work\\'             | '/c/Users/me/work/'
        'D:/data/reads.fq'                  | '/d/data/reads.fq'
        'C:\\'                              | '/c'
        'C:'                                | '/c'
        'C:\\My Data\\file (1).txt'         | '/c/My Data/file (1).txt'
        '\\\\?\\C:\\very\\long\\path'       | '/c/very/long/path'
        'C:file.txt'                        | 'C:file.txt'
    }

    def 'should map a local path object' () {
        expect:
        PosixPaths.toPosix((Path)null) == null
        PosixPaths.toPosix(Path.of('/some/dir')) == '/some/dir'
    }

    def 'should detect windows paths' () {
        expect:
        PosixPaths.isWindowsPath('C:\\data')
        PosixPaths.isWindowsPath('z:/data')
        PosixPaths.isWindowsPath('\\\\?\\C:\\data')
        !PosixPaths.isWindowsPath('/c/data')
        !PosixPaths.isWindowsPath('s3://bucket/data')
        !PosixPaths.isWindowsPath('data')
        !PosixPaths.isWindowsPath(null)
    }

    def 'should reject UNC paths' () {
        expect:
        PosixPaths.isUncPath('\\\\server\\share\\data')
        !PosixPaths.isUncPath('\\\\?\\C:\\data')
        !PosixPaths.isUncPath('/data')

        when:
        PosixPaths.toPosix('\\\\server\\share\\data')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('UNC paths are not supported')
    }

}
