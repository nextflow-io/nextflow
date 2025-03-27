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

package nextflow.data.cid.fs

import nextflow.data.cid.DefaultCidStore
import spock.lang.Shared

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
 * CID File system provider tests
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CidFileSystemProviderTest extends Specification {

    @Shared def wdir = Files.createTempDirectory('wdir')
    @Shared def meta = wdir.resolve('.meta')
    @Shared def data = wdir.resolve('work')

    def setupSpec(){
        meta.mkdirs()
        data.mkdirs()
    }

    def cleanupSpec(){
        wdir.deleteDir()
    }

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
        def config = [store:[location:data.toString()]]
        def cid = CidPath.asUri('cid://12345')
        when:
        def fs = provider.newFileSystem(cid, config) as CidFileSystem
        then:
        (fs.cidStore as DefaultCidStore).location == data
    }

    def 'should get a file system' () {
        given:
        def provider = new CidFileSystemProvider()
        def config = [store:[location: data.toString()]]
        def uri = CidPath.asUri('cid://12345')
        when:
        provider.getFileSystem(uri)
        then:
        thrown(FileSystemNotFoundException)

        when:
        provider.newFileSystem(uri, config) as CidFileSystem
        and:
        def fs = provider.getFileSystem(uri) as CidFileSystem
        then:
        (fs.cidStore as DefaultCidStore).location == data
    }

    def 'should get or create a file system' () {
        given:
        def config = [workflow:[data:[store:[location: data.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def uri = CidPath.asUri('cid://12345')
        def provider = new CidFileSystemProvider()
        
        when:
        def fs = provider.getFileSystemOrCreate(uri) as CidFileSystem
        then:
        (fs.cidStore as DefaultCidStore).location == data

        when:
        def fs2 = provider.getFileSystemOrCreate(uri) as CidFileSystem
        then:
        fs2.is(fs)
    }

    def 'should create new byte channel' () {
        given:
        def config = [workflow:[data:[store:[location:wdir.toString()]]]]
        def outputMeta = meta.resolve("12345/output.txt")
        def output = data.resolve("output.txt")
        output.text = "Hello, World!"
        outputMeta.mkdirs()
        outputMeta.resolve(".data.json").text = '{"type":"WorkflowOutput","path":"'+output.toString()+'"}'

        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath(CidPath.asUri('cid://12345/output.txt'))
        def opts = Set.of(StandardOpenOption.READ)
        when:
        def channel = provider.newByteChannel(cid, opts)
        and:
        def buffer = ByteBuffer.allocate(1000);
        def read = channel.read(buffer)
        channel.close()
        def bytes = new byte[read]
        buffer.get(0,bytes)
        then:
        bytes == "Hello, World!".getBytes()

        cleanup:
        outputMeta.deleteDir()
        output.delete()
    }

    def 'should read cid' () {
        given:
        def config = [workflow:[data:[store:[location:wdir.toString()]]]]
        def outputMeta = meta.resolve("12345/output.txt")
        def output = data.resolve("output.txt")
        output.text = "Hello, World!"
        outputMeta.mkdirs()
        outputMeta.resolve(".data.json").text = '{"type":"WorkflowOutput","path":"'+output.toString()+'"}'

        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath(CidPath.asUri('cid://12345/output.txt'))
        def opts = Set.of(StandardOpenOption.READ)

        expect:
        cid.text == "Hello, World!"

        cleanup:
        outputMeta.deleteDir()
        output.delete()
    }

    def 'should not create a directory' () {
        given:
        def config = [workflow:[data:[store:[location:'test']]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath(CidPath.asUri('cid://12345'))

        when:
        provider.createDirectory(cid)
        then:
        thrown(UnsupportedOperationException)

    }

    def 'should create directory stream' () {
        given:
        def output1 = data.resolve('path')
        output1.mkdir()
        output1.resolve('file1.txt').text = 'file1'
        output1.resolve('file2.txt').text = 'file2'
        output1.resolve('file3.txt').text = 'file3'
        meta.resolve('12345/output1').mkdirs()
        meta.resolve('12345/output2').mkdirs()
        meta.resolve('12345/.data.json').text = '{"type":"TaskRun"}'
        meta.resolve('12345/output1/.data.json').text = '{"type":"TaskOutput", "path": "' + output1.toString() + '"}'

        and:
        def config = [workflow:[data:[store:[location:wdir.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath(CidPath.asUri('cid://12345/output1'))
        def cid2 = provider.getPath(CidPath.asUri('cid://12345'))

        expect:
        Files.exists(cid)
        Files.exists(cid.resolve('file1.txt'))
        Files.exists(cid.resolve('file2.txt'))
        Files.exists(cid.resolve('file3.txt'))

        when:
        provider.newDirectoryStream(cid2, (p) -> true)
        then:
        thrown(FileNotFoundException)

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
        meta.resolve('12345').deleteDir()
        output1.deleteDir()

    }

    def 'should not delete a file' () {
        given:
        def config = [workflow:[data:[store:[location:'test']]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid = provider.getPath(CidPath.asUri('cid://12345'))
        
        when:
        provider.delete(cid)
        then:
        thrown(UnsupportedOperationException)

    }

    def 'should not copy a file' () {
        given:
        def config = [workflow:[data:[store:[location:'test']]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath(CidPath.asUri('cid://12345/abc'))
        def cid2 = provider.getPath(CidPath.asUri('cid://54321/foo'))

        when:
        provider.copy(cid1, cid2)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should not move a file' () {
        given:
        def config = [workflow:[data:[store:[location:'test']]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath(CidPath.asUri('cid://12345/abc'))
        def cid2 = provider.getPath(CidPath.asUri('cid://54321/foo'))

        when:
        provider.move(cid1, cid2)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should check is same file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [workflow:[data:[store:[location:folder.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath(CidPath.asUri('cid://12345/abc'))
        def cid2 = provider.getPath(CidPath.asUri('cid://54321/foo'))
        def cid3 = provider.getPath(CidPath.asUri('cid://54321/foo'))

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
        def config = [workflow:[data:[store:[location:wdir.toString()]]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def output = folder.resolve('path')
        output.mkdir()
        output.resolve('abc').text = 'file1'
        output.resolve('.foo').text = 'file2'
        meta.resolve('12345/output').mkdirs()
        meta.resolve('12345/output/.data.json').text = '{"type":"TaskOutput", "path": "' + output.toString() + '"}'
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath(CidPath.asUri('cid://12345/output/abc'))
        def cid2 = provider.getPath(CidPath.asUri('cid://12345/output/.foo'))

        expect:
        !provider.isHidden(cid1)
        provider.isHidden(cid2)

        cleanup:
        folder?.deleteDir()
    }

    def 'should read file attributes' () {
        given:
        def config = [workflow:[data:[store:[location:wdir.toString()]]]]
        def file = data.resolve('abc')
        file.text = 'Hello'
        meta.resolve('12345/abc').mkdirs()
        meta.resolve('12345/abc/.data.json').text = '{"type":"TaskOutput", "path": "' + file.toString() + '"}'
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new CidFileSystemProvider()
        def cid1 = provider.getPath(CidPath.asUri('cid://12345/abc'))

        when:
        def attr1 = provider.readAttributes(cid1, BasicFileAttributes)
        def real1= Files.readAttributes(file,BasicFileAttributes)
        then:
        !attr1.directory
        attr1.isRegularFile()
        attr1.size() == real1.size()
        attr1.creationTime() == real1.creationTime()
        attr1.lastModifiedTime() == real1.lastModifiedTime()
        attr1.lastAccessTime() == real1.lastAccessTime()

        cleanup:
        file?.delete()
        meta.resolve('12345').deleteDir()
    }

}

