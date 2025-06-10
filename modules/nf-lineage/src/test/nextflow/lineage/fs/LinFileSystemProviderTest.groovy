/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.fs

import java.nio.ByteBuffer
import java.nio.channels.NonWritableChannelException
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes

import nextflow.Global
import nextflow.Session
import nextflow.lineage.DefaultLinStore
import spock.lang.Shared
import spock.lang.Specification
/**
 * LID File system provider tests
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinFileSystemProviderTest extends Specification {

    @Shared def wdir = Files.createTempDirectory('wdir')
    @Shared def data = Files.createTempDirectory('work')


    def cleanupSpec(){
        wdir.deleteDir()
        data.deleteDir()
    }

    def 'should return lid scheme' () {
        given:
        def provider = new LinFileSystemProvider()
        expect:
        provider.getScheme() == 'lid'
    }

    def 'should get lid path' () {
        given:
        def lid = Mock(LinPath)
        and:
        def provider = new LinFileSystemProvider()
        expect:
        provider.toLinPath(lid) == lid

        when:
        provider.toLinPath(Path.of('foo'))
        then:
        thrown(ProviderMismatchException)
    }

    def 'should create new file system' () {
        given:
        def provider = new LinFileSystemProvider()
        def config = [store:[location:data.toString()]]
        def lid = LinPath.asUri('lid://12345')
        when:
        def fs = provider.newFileSystem(lid, config) as LinFileSystem
        then:
        (fs.store as DefaultLinStore).location == data
    }

    def 'should get a file system' () {
        given:
        def provider = new LinFileSystemProvider()
        def config = [store:[location: data.toString()]]
        def uri = LinPath.asUri('lid://12345')
        when:
        provider.getFileSystem(uri)
        then:
        thrown(FileSystemNotFoundException)

        when:
        provider.newFileSystem(uri, config) as LinFileSystem
        and:
        def fs = provider.getFileSystem(uri) as LinFileSystem
        then:
        (fs.store as DefaultLinStore).location == data
    }

    def 'should get or create a file system' () {
        given:
        def config = [lineage:[store:[location: data.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def uri = LinPath.asUri('lid://12345')
        def provider = new LinFileSystemProvider()
        
        when:
        def fs = provider.getFileSystemOrCreate(uri) as LinFileSystem
        then:
        (fs.store as DefaultLinStore).location == data

        when:
        def fs2 = provider.getFileSystemOrCreate(uri) as LinFileSystem
        then:
        fs2.is(fs)
    }

    def 'should create new byte channel' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        def outputMeta = wdir.resolve("12345/output.txt")
        def output = data.resolve("output.txt")
        output.text = "Hello, World!"
        outputMeta.mkdirs()
        outputMeta.resolve(".data.json").text = '{"version":"lineage/v1beta1","kind":"FileOutput","path":"'+output.toString()+'"}'

        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345/output.txt'))
        def opts = Set.of(StandardOpenOption.READ)
        when:
        def channel = provider.newByteChannel(lid, opts)
        then:
        channel.isOpen()
        channel.position() == 0
        channel.size() == "Hello, World!".getBytes().size()
        when:
        channel.truncate(25)
        then:
        thrown(NonWritableChannelException)

        when:
        def buffer = ByteBuffer.allocate(1000);
        def read = channel.read(buffer)
        def bytes = new byte[read]
        buffer.get(0,bytes)
        then:
        bytes == "Hello, World!".getBytes()
        when:
        channel.position(2)
        then:
        channel.position() == 2

        when:
        channel.write(buffer)
        then:
        thrown(NonWritableChannelException)

        when:
        provider.newByteChannel(lid, Set.of(StandardOpenOption.WRITE))
        then:
        thrown(UnsupportedOperationException)

        when:
        provider.newByteChannel(lid, Set.of(StandardOpenOption.APPEND))
        then:
        thrown(UnsupportedOperationException)

        cleanup:
        channel.close()
        outputMeta.deleteDir()
        output.delete()
    }

    def 'should create new byte channel for LinMetadata' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        def outputMeta = wdir.resolve("12345")
        outputMeta.mkdirs()
        outputMeta.resolve(".data.json").text = '{"version":"lineage/v1beta1","kind":"WorkflowRun","sessionId":"session","name":"run_name","params":[{"type":"String","name":"param1","value":"value1"}]}'

        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345#name'))

        when:
        def channel = provider.newByteChannel(lid, Set.of(StandardOpenOption.READ))
        then:
        channel.isOpen()
        channel.position() == 0
        channel.size() == '"run_name"'.getBytes().size()

        when:
        channel.truncate(25)
        then:
        thrown(NonWritableChannelException)

        when:
        def buffer = ByteBuffer.allocate(1000);
        def read = channel.read(buffer)
        def bytes = new byte[read]
        buffer.get(0,bytes)
        then:
        bytes =='"run_name"'.getBytes()

        when:
        channel.position(2)
        then:
        channel.position() == 2

        when:
        channel.write(buffer)
        then:
        thrown(NonWritableChannelException)

        when:
        provider.newByteChannel(lid, Set.of(StandardOpenOption.WRITE))
        then:
        thrown(UnsupportedOperationException)

        when:
        provider.newByteChannel(lid, Set.of(StandardOpenOption.APPEND))
        then:
        thrown(UnsupportedOperationException)

        cleanup:
        channel.close()
        outputMeta.deleteDir()
    }

    def 'should read lid' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        def outputMeta = wdir.resolve("12345/output.txt")
        def output = data.resolve("output.txt")
        output.text = "Hello, World!"
        outputMeta.mkdirs()
        outputMeta.resolve(".data.json").text = '{"version":"lineage/v1beta1","kind":"FileOutput","path":"'+output.toString()+'"}'

        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345/output.txt'))
        def opts = Set.of(StandardOpenOption.READ)

        expect:
        lid.text == "Hello, World!"

        cleanup:
        outputMeta.deleteDir()
        output.delete()
    }

    def 'should not create a directory' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345'))

        when:
        provider.createDirectory(lid)
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
        wdir.resolve('12345/output1').mkdirs()
        wdir.resolve('12345/output2').mkdirs()
        wdir.resolve('12345/.data.json').text = '{"version":"lineage/v1beta1","kind":"TaskRun"}'
        wdir.resolve('12345/output1/.data.json').text = '{"version":"lineage/v1beta1","kind":"FileOutput", "path": "' + output1.toString() + '"}'

        and:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345/output1'))
        def lid2 = provider.getPath(LinPath.asUri('lid://12345'))
        def lid3 = provider.getPath(LinPath.asUri('lid://678'))

        expect:
        Files.exists(lid)
        Files.exists(lid.resolve('file1.txt'))
        Files.exists(lid.resolve('file2.txt'))
        Files.exists(lid.resolve('file3.txt'))

        when:
        provider.newDirectoryStream(lid3, (p) -> true)
        then:
        thrown(FileNotFoundException)

        when:
        def stream = provider.newDirectoryStream(lid, (p) -> true)
        and:
        def result = stream.toList()
        then:
        result.toSet() == [
            lid.resolve('file1.txt'),
            lid.resolve('file2.txt'),
            lid.resolve('file3.txt')
        ] as Set

        when:
        stream = provider.newDirectoryStream(lid2, (p) -> true)
        and:
        result = stream.toList()
        then:
        result.size() == 1
        result[0] ==  lid2.resolve('output1')

        cleanup:
        wdir.resolve('12345').deleteDir()
        output1.deleteDir()

    }

    def 'should not delete a file' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid = provider.getPath(LinPath.asUri('lid://12345'))
        
        when:
        provider.delete(lid)
        then:
        thrown(UnsupportedOperationException)

    }

    def 'should not copy a file' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/abc'))
        def lid2 = provider.getPath(LinPath.asUri('lid://54321/foo'))

        when:
        provider.copy(lid1, lid2)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should not move a file' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/abc'))
        def lid2 = provider.getPath(LinPath.asUri('lid://54321/foo'))

        when:
        provider.move(lid1, lid2)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should check is same file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage:[store:[location:folder.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/abc'))
        def lid2 = provider.getPath(LinPath.asUri('lid://54321/foo'))
        def lid3 = provider.getPath(LinPath.asUri('lid://54321/foo'))

        expect:
        !provider.isSameFile(lid1, lid2)
        !provider.isSameFile(lid1, lid3)
        and:
        provider.isSameFile(lid2, lid3)

        cleanup:
        folder?.deleteDir()
    }

    def 'should check is hidden file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def output = folder.resolve('path')
        output.mkdir()
        output.resolve('abc').text = 'file1'
        output.resolve('.foo').text = 'file2'
        wdir.resolve('12345/output').mkdirs()
        wdir.resolve('12345/output/.data.json').text = '{"version":"lineage/v1beta1","kind":"FileOutput", "path": "' + output.toString() + '"}'
        and:
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/output/abc'))
        def lid2 = provider.getPath(LinPath.asUri('lid://12345/output/.foo'))

        expect:
        !provider.isHidden(lid1)
        provider.isHidden(lid2)

        cleanup:
        folder?.deleteDir()
    }

    def 'should read file attributes' () {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        def file = data.resolve('abc')
        file.text = 'Hello'
        wdir.resolve('12345/abc').mkdirs()
        wdir.resolve('12345/abc/.data.json').text = '{"version":"lineage/v1beta1","kind":"FileOutput", "path":"' + file.toString() + '"}'
        and:
        Global.session = Mock(Session) { getConfig()>>config }
        and:
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/abc'))

        when:
        def attr1 = provider.readAttributes(lid1, BasicFileAttributes)
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
        wdir.resolve('12345').deleteDir()
    }

    def 'should throw exception in unsupported methods'() {
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        def provider = new LinFileSystemProvider()

        when:
        provider.newOutputStream(null)
        then:
        thrown(UnsupportedOperationException)

        when:
        provider.getFileStore(null)
        then:
        thrown(UnsupportedOperationException)

        when:
        provider.readAttributes(null, "attrib")
        then:
        thrown(UnsupportedOperationException)

        when:
        provider.setAttribute(null, "attrib", null)
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should throw exception when checking access mode'(){
        given:
        def config = [lineage:[store:[location:wdir.toString()]]]
        Global.session = Mock(Session) { getConfig()>>config }
        def provider = new LinFileSystemProvider()
        def lid1 = provider.getPath(LinPath.asUri('lid://12345/abc'))

        when:
        provider.checkAccess(lid1, AccessMode.WRITE)
        then:
        def ex1 = thrown(AccessDeniedException)
        ex1.message == "Write mode not supported"

        when:
        provider.checkAccess(lid1, AccessMode.EXECUTE)
        then:
        def ex2 = thrown(AccessDeniedException)
        ex2.message == "Execute mode not supported"
    }
}

