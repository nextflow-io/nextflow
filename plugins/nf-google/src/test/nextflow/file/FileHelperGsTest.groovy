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

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
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
                CloudStorageFileSystem.forBucket('foo').getPath('')

        and:
        FileHelper.asPath('gs://foo/this/and/that.txt') ==
                CloudStorageFileSystem.forBucket('foo').getPath('/this/and/that.txt')

        and:
        FileHelper.asPath('gs://foo/b a r.txt') ==
                CloudStorageFileSystem.forBucket('foo').getPath('/b a r.txt')

        and:
        FileHelper.asPath('gs://f o o/bar.txt') ==
                CloudStorageFileSystem.forBucket('f o o').getPath('/bar.txt')

        and:
        FileHelper.asPath('gs://f_o_o/bar.txt') ==
                CloudStorageFileSystem.forBucket('f_o_o').getPath('/bar.txt')
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
        def foo = CloudStorageFileSystem.forBucket('nf-bucket').getPath('foo.txt')
        def bar = CloudStorageFileSystem.forBucket('nf-bucket').getPath('bar.txt')
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
        FileHelper.toCanonicalPath(VALUE) == (EXPECTED ? FileHelper.asPath(EXPECTED) : null)

        cleanup:
        SysEnv.pop()
        Global.session = null

        where:
        VALUE                       | EXPECTED
        null                        | null
        'file.txt'                  | 'gs://host.com/work/file.txt'
        Path.of('file.txt')         | 'gs://host.com/work/file.txt'
        and:
        './file.txt'                | 'gs://host.com/work/file.txt'
        '.'                         | 'gs://host.com/work'
        './'                        | 'gs://host.com/work'
        '../file.txt'               | 'gs://host.com/file.txt'
        and:
        '/file.txt'                 | '/file.txt'
        Path.of('/file.txt')        | '/file.txt'
    }

    def 'should convert to a canonical path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [google:[project:'foo', region:'x']] }

        expect:
        FileHelper.toCanonicalPath(VALUE).toUri() == EXPECTED

        where:
        VALUE                       | EXPECTED
        'gs://foo/some/file.txt'    | new URI('gs://foo/some/file.txt')
        'gs://foo/some///file.txt'  | new URI('gs://foo/some/file.txt')
    }

    @Unroll
    def 'should remove consecutive slashes in the path' () {
        given:
        Global.session = Mock(Session) { getConfig() >> [google:[project:'foo', region:'x']] }

        expect:
        FileHelper.asPath(STR).toUri() == EXPECTED
        where:
        STR                         | EXPECTED
        'gs://foo//this/that'       | new URI('gs://foo/this/that')
        'gs://foo//this///that'     | new URI('gs://foo/this/that')
    }

}
