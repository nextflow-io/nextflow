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
 */

package nextflow.file

import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.google.nio.GsFileSystem
import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileHelperGsTest extends Specification {

    def 'should parse google storage path' () {

        given:
        Global.session = Mock(Session) {
            getConfig() >> [google:[project:'foo', region:'x']]
        }

        expect:
        FileHelper.asPath('file.txt') ==
                Paths.get('file.txt')
        and:
        FileHelper.asPath('gs://foo') ==
                GsFileSystem.forBucket('foo').getPath('')

        and:
        FileHelper.asPath('gs://foo/this/and/that.txt') ==
                GsFileSystem.forBucket('foo').getPath('/this/and/that.txt')

        and:
        FileHelper.asPath('gs://foo/b a r.txt') ==
                GsFileSystem.forBucket('foo').getPath('/b a r.txt')

        and:
        FileHelper.asPath('gs://f o o/bar.txt') ==
                GsFileSystem.forBucket('f o o').getPath('/bar.txt')

        and:
        FileHelper.asPath('gs://f_o_o/bar.txt') ==
                GsFileSystem.forBucket('f_o_o').getPath('/bar.txt')
    }


    def 'should strip ending slash' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [google:[project:'foo', region:'x']] }
        def nxFolder = Paths.get('/my-bucket/foo')
        def nxNested = Paths.get('/my-bucket/foo/bar/')
        and:
        def gsFolder = 'gs://my-bucket/foo' as Path
        def gsNested = 'gs://my-bucket/foo/bar/' as Path

        expect:
        nxFolder.relativize(nxNested).toString() == 'bar'
        gsFolder.relativize(gsNested).toString() == 'bar/'      // <-- gs adds a slash that mess-up things
        and:
        FileHelper.relativize0(nxFolder,nxNested).toString() == 'bar'
        FileHelper.relativize0(gsFolder,gsNested).toString() == 'bar'

    }

    @Ignore
    def 'should throw FileAlreadyExistsException'() {
        given:
        def foo = GsFileSystem.forBucket('nf-bucket').getPath('foo.txt')
        def bar = GsFileSystem.forBucket('nf-bucket').getPath('bar.txt')
        and:
        if( !Files.exists(foo) ) Files.createFile(foo)
        if( !Files.exists(bar) ) Files.createFile(bar)
        when:
        Files.copy(foo, bar)
        then:
        thrown(FileAlreadyExistsException)
    }

    @Unroll
    def 'should convert to canonical path with base' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [google:[project:'foo', region:'x']] }
        and:
        SysEnv.push(NXF_FILE_ROOT: 'gs://host.com/work')

        expect:
        FileHelper.toCanonicalPath(VALUE) == EXPECTED

        cleanup:
        SysEnv.pop()
        Global.session = null

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | FileHelper.asPath('gs://host.com/work/file.txt')
        Path.of('file.txt')         | FileHelper.asPath('gs://host.com/work/file.txt')
        and:
        './file.txt'                | FileHelper.asPath('gs://host.com/work/file.txt')
        '.'                         | FileHelper.asPath('gs://host.com/work')
        './'                        | FileHelper.asPath('gs://host.com/work')
        '../file.txt'               | FileHelper.asPath('gs://host.com/file.txt')
        and:
        '/file.txt'                 | Path.of('/file.txt')
        Path.of('/file.txt')        | Path.of('/file.txt')
    }
}
