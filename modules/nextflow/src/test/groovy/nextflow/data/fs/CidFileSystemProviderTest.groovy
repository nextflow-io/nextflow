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

package nextflow.data.fs

import java.nio.ByteBuffer
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes

import nextflow.Global
import nextflow.Session
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidFileSystemProviderTest extends Specification {

    def 'should return cid scheme' () {
        given:
        def provider = new CidFileSystemProvider()
        expect:
        provider.getScheme() == 'cid'
    }

    def 'should get cid path' () {
        given:
        def cid = Mock(CidPath)
        and:
        def provider = new CidFileSystemProvider()
        expect:
        provider.toCidPath(cid) == cid

        when:
        provider.toCidPath(Path.of('foo'))
        then:
        thrown(ProviderMismatchException)
    }

    def 'should create new file system' () {
        given:
        def provider = new CidFileSystemProvider()
        def config = [store:[location:'/data']]
        def cid = CidPath.asUri('cid://12345')
        when:
        def fs = provider.newFileSystem(cid, config) as CidFileSystem
        then:
        fs.basePath == Path.of('/data')
    }

    def 'should get a file system' () {
        given:
        def provider = new CidFileSystemProvider()
        def config = [store:[location:'/data']]
        def uri = CidPath.asUri('cid://12345')
        when:
        provider.getFileSystem(uri)
        then:
        thrown(FileSystemNotFoundException)

        when:
        provider.newFileSystem(uri, config) as CidFileSystem
        and:
        def result = provider.getFileSystem(uri) as CidFileSystem
        then:
        result.basePath == Path.of('/data')
    }

    def 'should get or create a file system' () {
        given:
        def config = [workflow:[data:[store:[location:'/this/that']]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def uri = CidPath.asUri('cid://12345')
        def provider = new CidFileSystemProvider()
        
        when:
        def fs = provider.getFileSystemOrCreate(uri) as CidFileSystem
        then:
        fs.basePath == Path.of('/this/that')

        when:
        def fs2 = provider.getFileSystemOrCreate(uri) as CidFileSystem
        then:
        fs2.is(fs)
    }

    def 'should get a path' () {
       given:
       def config = [workflow:[data:[store:[location:'/data']]]]
       Global.session = Mock(Session) { getConfig()>>config }
       and:
       def provider = new CidFileSystemProvider()
       def uri1 = CidPath.asUri('cid://12345')
       def uri2 = CidPath.asUri('cid://12345/foo/bar')

        when:
        def cid1 = provider.getPath(uri1)
        then:
        cid1.toRealPath() == Path.of('/data/12345')

        when:
        def cid2 = provider.getPath(uri2)
        then:
        cid2.toRealPath() == Path.of('/data/12345/foo/bar')
    }

    def 'should create new byte channel' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath('cid://12345')
        def opts = Set.of(StandardOpenOption.WRITE, StandardOpenOption.CREATE)

        when:
        def channel = provider.newByteChannel(cid, opts)
        and:
        def buffer = ByteBuffer.wrap("Hello, World!".getBytes());
        channel.write(buffer)
        channel.close()
        then:
        cid.text == "Hello, World!"

        cleanup:
        folder?.deleteDir()
    }

    def 'should create a directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath('cid://12345')

        when:
        provider.createDirectory(cid)
        then:
        Files.isDirectory(cid)

        cleanup:
        folder?.deleteDir()
    }

    def 'should create many dirs' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345')
        def cid2 = provider.getPath('cid://54321/foo/bar')

        when:
        Files.createDirectories(cid1)
        then:
        noExceptionThrown()

        when:
        Files.createDirectories(cid2)
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()
    }

    def 'should create directory stream' () {
        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('12345').mkdir()
        folder.resolve('12345/file1.txt').text = 'file1'
        folder.resolve('12345/file2.txt').text = 'file2'
        folder.resolve('12345/file3.txt').text = 'file3'
        and:
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath('cid://12345')

        expect:
        Files.exists(cid)
        Files.exists(cid.resolve('file1.txt'))
        Files.exists(cid.resolve('file2.txt'))
        Files.exists(cid.resolve('file3.txt'))

        when:
        def stream = provider.newDirectoryStream(cid, (p) -> true)
        and:
        def result = stream.toList()
        then:
        result.toSet() == [
            cid.resolve('file1.txt'),
            cid.resolve('file2.txt'),
            cid.resolve('file3.txt')
        ] as Set

        cleanup:
        folder?.deleteDir()
    }

    def 'should delete a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath('cid://12345')

        when:
        cid.text = 'new file'
        then:
        Files.exists(cid)
        
        when:
        provider.delete(cid)
        then:
        !Files.exists(cid)

        cleanup:
        folder?.deleteDir()
    }

    def 'should copy a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345/abc')
        def cid2 = provider.getPath('cid://54321/foo')
        and:
        Files.createDirectories(cid1.parent)
        Files.createDirectories(cid2.parent)
        cid1.text = 'some data'

        when:
        provider.copy(cid1, cid2)
        then:
        Files.exists(cid1)
        Files.exists(cid2)
        cid2.text == cid1.text

        cleanup:
        folder?.deleteDir()
    }

    def 'should move a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345/abc')
        def cid2 = provider.getPath('cid://54321/foo')
        and:
        Files.createDirectories(cid1.parent)
        Files.createDirectories(cid2.parent)
        cid1.text = 'some data'

        when:
        provider.move(cid1, cid2)
        then:
        !Files.exists(cid1)
        Files.exists(cid2)
        cid2.text == 'some data'

        cleanup:
        folder?.deleteDir()
    }

    def 'should check is same file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345/abc')
        def cid2 = provider.getPath('cid://54321/foo')
        def cid3 = provider.getPath('cid://54321/foo')

        expect:
        !provider.isSameFile(cid1, cid2)
        !provider.isSameFile(cid1, cid3)
        and:
        provider.isSameFile(cid2, cid3)

        cleanup:
        folder?.deleteDir()
    }

    def 'should check is hidden file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345/abc')
        def cid2 = provider.getPath('cid://54321/.foo')

        expect:
        !provider.isHidden(cid1)
        provider.isHidden(cid2)

        cleanup:
        folder?.deleteDir()
    }

    def 'should read file attributes' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath('cid://12345/abc')
        and:
        Files.createDirectories(cid1.parent)
        cid1.text = 'Hello'

        when:
        def attr1 = provider.readAttributes(cid1, BasicFileAttributes)
        def real1= Files.readAttributes(cid1.toRealPath(),BasicFileAttributes)
        then:
        attr1.fileKey() == cid1.fileId
        and:
        !attr1.directory
        attr1.isRegularFile()
        attr1.size() == real1.size()
        attr1.creationTime() == real1.creationTime()
        attr1.lastModifiedTime() == real1.lastModifiedTime()
        attr1.lastAccessTime() == real1.lastAccessTime()

        when:
        def cid2 = cid1.parent as CidPath
        def attr2 = provider.readAttributes(cid2, BasicFileAttributes)
        def real2= Files.readAttributes(cid2.toRealPath(),BasicFileAttributes)
        then:
        attr2.fileKey() == cid2.fileId
        and:
        attr2.directory
        !attr2.isRegularFile()
        attr2.size() == real2.size()
        attr2.creationTime() == real2.creationTime()
        attr2.lastModifiedTime() == real2.lastModifiedTime()
        attr2.lastAccessTime() == real2.lastAccessTime()

        cleanup:
        folder?.deleteDir()
    }

}

