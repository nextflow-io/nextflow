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

package nextflow.cloud.aws.nio

import java.nio.charset.Charset
import java.nio.file.DirectoryNotEmptyException
import java.nio.file.FileAlreadyExistsException
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.SimpleFileVisitor
import java.nio.file.StandardCopyOption
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes

import com.amazonaws.services.s3.AmazonS3
import com.amazonaws.services.s3.model.AmazonS3Exception
import com.amazonaws.services.s3.model.Tag
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.file.CopyMoveHelper
import nextflow.file.FileHelper
import nextflow.trace.TraceHelper
import spock.lang.Ignore
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Timeout
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Timeout(60)
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AWS_S3FS_ACCESS_KEY') && System.getenv('AWS_S3FS_SECRET_KEY')})
class AwsS3NioTest extends Specification implements AwsS3BaseSpec {

    @Shared
    static AmazonS3 s3Client0

    AmazonS3 getS3Client() { s3Client0 }

    static {
        def fs = (S3FileSystem)FileHelper.getOrCreateFileSystemFor(URI.create("s3:///"), config0())
        s3Client0 = fs.client.getClient()
    }

    static private Map config0() {
        def accessKey = System.getenv('AWS_S3FS_ACCESS_KEY')
        def secretKey = System.getenv('AWS_S3FS_SECRET_KEY')
        return [aws: [access_key: accessKey, secret_key: secretKey]]
    }

    def setup() {
        def cfg = config0()
        Global.config = cfg
        Global.session = Mock(Session) { getConfig()>>cfg }
    }

    def 'should create a blob' () {
        given:
        def bucket = createBucket()
        def path = s3path("s3://$bucket/file-name.txt")

        when:
        Files.createFile(path)
        then:
        existsPath("$bucket/file-name.txt")

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    def 'should write a file' () {
        given:
        def TEXT = "Hello world!"
        and:
        def bucket = createBucket()
        def path = s3path("s3://$bucket/file-name.txt")

        when:
        Files.write(path, TEXT.bytes)
        then:
        existsPath("$bucket/file-name.txt")
        readObject(path) == TEXT

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    def 'should read a file' () {
        given:
        def TEXT = "Hello world!"
        and:
        def bucket = createBucket()
        def path = s3path("s3://$bucket/file-name.txt")
        when:
        createObject("$bucket/file-name.txt", TEXT)
        then:
        new String(Files.readAllBytes(path)) == TEXT
        Files.readAllLines(path, Charset.forName('UTF-8')).get(0) == TEXT

        cleanup:
        if( bucket ) deleteBucket(bucket)
    }

    def 'should read file attributes' () {
        given:
        final start = System.currentTimeMillis()
        final TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def objectKey = 'data/alpha.txt'
        def filePath = "$bucketName/$objectKey"
        createObject(filePath, TEXT)
        and:

        //
        // -- readAttributes
        //
        def path = s3path("s3://$filePath")
        def attrs = Files.readAttributes(path, BasicFileAttributes)
        then:
        attrs.isRegularFile()
        !attrs.isDirectory()
        attrs.size() == 12
        !attrs.isSymbolicLink()
        !attrs.isOther()
        attrs.fileKey() == objectKey
        attrs.lastAccessTime().toMillis()-start < 5_000
        attrs.lastModifiedTime().toMillis()-start < 5_000
        attrs.creationTime().toMillis()-start < 5_000

        //
        // -- getLastModifiedTime
        //
        when:
        def time = Files.getLastModifiedTime(path)
        then:
        time == attrs.lastModifiedTime()

        //
        // -- getFileAttributeView
        //
// FIXME not supported
//        when:
//        def view = Files.getFileAttributeView(path, BasicFileAttributeView)
//        then:
//        view.readAttributes() == attrs

        //
        // -- readAttributes for a directory
        //
        when:
        attrs = Files.readAttributes(path.getParent(), BasicFileAttributes)
        then:
        !attrs.isRegularFile()
        attrs.isDirectory()
        attrs.size() == 0
        !attrs.isSymbolicLink()
        !attrs.isOther()
        attrs.fileKey() == "data/"
        attrs.lastAccessTime() .toMillis()-start < 5_000
        attrs.lastModifiedTime() .toMillis()-start < 5_000
        attrs.creationTime() .toMillis()-start < 5_000

        //
        // -- readAttributes for a bucket
        //
        when:
        attrs = Files.readAttributes(s3path("s3://$bucketName"), BasicFileAttributes)
        then:
        !attrs.isRegularFile()
        attrs.isDirectory()
        attrs.size() == 0
        !attrs.isSymbolicLink()
        !attrs.isOther()
        attrs.fileKey() == "/"
        attrs.creationTime() == null
        attrs.lastAccessTime() == null
        attrs.lastModifiedTime() == null

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }


    def 'should copy a stream to bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/data/file.txt")
        
        and:
        def stream = new ByteArrayInputStream(new String(TEXT).bytes)
        Files.copy(stream, target)
        then:
        existsPath(target)
        readObject(target) == TEXT

        when:
        stream = new ByteArrayInputStream(new String(TEXT).bytes)
        Files.copy(stream, target, StandardCopyOption.REPLACE_EXISTING)
        then:
        existsPath(target)
        readObject(target) == TEXT

        when:
        stream = new ByteArrayInputStream(new String(TEXT).bytes)
        Files.copy(stream, target)
        then:
        thrown(FileAlreadyExistsException)

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }


    def 'copy local file to a bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/data/file.txt")
        def source = Files.createTempFile('test','nf')
        source.text = TEXT

        and:
        Files.copy(source, target)
        then:
        readObject(target) == TEXT

        cleanup:
        if( source ) Files.delete(source)
        if( bucketName ) deleteBucket(bucketName)
    }

    def 'copy a remote file to a bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/target/file.txt")

        and:
        final objectName = "$bucketName/source/file.txt"
        final source = s3path("s3://$objectName")
        createObject(objectName, TEXT)

        and:
        Files.copy(source, target)
        then:
        existsPath(source)
        existsPath(target)
        readObject(target) == TEXT

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }

    @Ignore // FIXME
    def 'move a remote file to a bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/target/file.txt")

        and:
        final objectName = "$bucketName/source/file.txt"
        final source = s3path("s3://$objectName")
        createObject(objectName, TEXT)

        and:
        Files.move(source, target)
        then:
        !existsPath(source)
        existsPath(target)
        readObject(target) == TEXT

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }

    def 'move local file to a bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/target/file.txt")

        and:
        final source = Files.createTempFile('foo',null); source.text = TEXT

        and:
        Files.move(source, target)
        then:
        !Files.exists(source)
        existsPath(target)
        readObject(target) == TEXT

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }

    def 'move a remote file to local' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = Files.createTempFile('foo',null)

        and:
        final objectName = "$bucketName/source/file.txt"
        final source = s3path("s3://$objectName")
        createObject(objectName, TEXT)

        and:
        Files.move(source, target, StandardCopyOption.REPLACE_EXISTING)
        then:
        !existsPath(source)
        Files.exists(target)
        target.text == TEXT

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
        if( target ) Files.deleteIfExists(target)
    }

    @Ignore //FIXME
    def 'should create a directory' () {

        given:
        def bucketName = getRndBucketName()
        def dir = s3path("s3://$bucketName")

        when:
        Files.createDirectory(dir)
        then:
        existsPath(dir)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create a directory tree' () {
        given:
        def bucketName = createBucket()
        def dir = s3path("s3://$bucketName/alpha/bravo/omega/")
        when:
        Files.createDirectories(dir)
        then:
        Files.exists(s3path("s3://$bucketName/alpha/"))
        Files.exists(s3path("s3://$bucketName/alpha/bravo/"))
        Files.exists(s3path("s3://$bucketName/alpha/bravo/omega/"))

        when:
        Files.createDirectories(dir)
        then:
        noExceptionThrown()

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create a file' () {
        given:
        def bucketName = createBucket()

        when:
        def path = s3path("s3://$bucketName/data/file.txt")
        Files.createFile(path)
        then:
        existsPath(path)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create temp file and directory' () {
        given:
        def bucketName = createBucket()
        def base = s3path("s3://$bucketName")

        when:
        def t1 = Files.createTempDirectory(base, 'test')
        then:
        Files.exists(t1)

        when:
        def t2 = Files.createTempFile(base, 'prefix', 'suffix')
        then:
        Files.exists(t2)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should delete a file' () {
        given:
        def bucketName = createBucket()
        def target = s3path("s3://$bucketName/data/file.txt")
        and:
        createObject(target.toString(), 'HELLO WORLD')

        when:
        Files.delete(target)
        sleep 100
        then:
        !existsPath(target)

        cleanup:
        deleteBucket(bucketName)
    }

    @Ignore // FIXME
    def 'should delete a bucket' () {
        given:
        final bucketName = createBucket()

        when:
        Files.delete(s3path("s3://$bucketName"))
        then:
        !existsPath(bucketName)

    }

    @Ignore // FIXME
    def 'should throw when deleting a not empty container' () {
        given:
        def bucketName = createBucket()
        and:
        createObject("$bucketName/this/that", 'HELLO')

        when:
        def path1 = s3path("s3://$bucketName")
        Files.delete(path1)
        then:
        thrown(DirectoryNotEmptyException)

        when:
        def path2 = s3path("s3://$bucketName/this")
        Files.delete(path2)
        then:
        thrown(DirectoryNotEmptyException)

        when:
        createObject("$bucketName/this", 'HELLO')
        Files.delete(Paths.get(path2))
        then:
        thrown(DirectoryNotEmptyException)

        cleanup:
        deleteBucket(bucketName)
    }

    @Ignore // FIXME
    def 'should throw a NoSuchFileException when deleting an object not existing' () {

        given:
        def bucketName = getRndBucketName()
        def path = s3path("s3://$bucketName/alpha/bravo")

        when:
        Files.delete(path)
        then:
        thrown(NoSuchFileException)

    }

    @Ignore //FIXME
    def 'should validate exists method' () {
        given:
        def bucketName = createBucket()
        and:
        def missingBucket = getRndBucketName()
        and:
        createObject("$bucketName/file.txt", 'HELLO')

        expect:
        Files.exists(s3path("s3://$bucketName"))
        Files.exists(s3path("s3://$bucketName/file.txt"))
        !Files.exists(s3path("s3://$bucketName/fooooo.txt"))
        !Files.exists(s3path("s3://$missingBucket"))

        cleanup:
        deleteBucket(bucketName)
    }


    def 'should check is it is a directory' () {
        given:
        def bucketName = createBucket()

        when:
        def path = s3path("s3://$bucketName")
        then:
        Files.isDirectory(path)
        !Files.isRegularFile(path)

        when:
        def file = path.resolve('this/and/that')
        createObject(file, 'Hello world')
        then:
        !Files.isDirectory(file)
        Files.isRegularFile(file)
        Files.isReadable(file)
        Files.isWritable(file)
        !Files.isExecutable(file)
        !Files.isSymbolicLink(file)

        expect:
        Files.isDirectory(file.parent)
        !Files.isRegularFile(file.parent)
        Files.isReadable(file)
        Files.isWritable(file)
        !Files.isExecutable(file)
        !Files.isSymbolicLink(file)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should check that is the same file' () {

        given:
        def file1 = s3path("s3://some/data/file.txt")
        def file2 = s3path("s3://some/data/file.txt")
        def file3 = s3path("s3://some/data/fooo.txt")

        expect:
        Files.isSameFile(file1, file2)
        !Files.isSameFile(file1, file3)

    }

    def 'should create a newBufferedReader' () {
        given:
        def bucketName = createBucket()
        and:
        def TEXT = randomText(50 * 1024)
        def path = s3path("s3://$bucketName/file.txt")
        createObject(path, TEXT)

        when:
        def reader = Files.newBufferedReader(path, Charset.forName('UTF-8'))
        then:
        reader.text == TEXT

        when:
        def unknown = s3path("s3://$bucketName/unknown.txt")
        Files.newBufferedReader(unknown, Charset.forName('UTF-8'))
        then:
        thrown(NoSuchFileException)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create a newBufferedWriter' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(50 * 1024)
        final path = s3path("s3://$bucketName/file.txt")

        when:
        def writer = Files.newBufferedWriter(path, Charset.forName('UTF-8'))
        TEXT.readLines().each { it -> writer.println(it) }
        writer.close()
        then:
        readObject(path) == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create a newInputStream' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(50 * 1024)
        final path = s3path("s3://$bucketName/file.txt")
        createObject(path, TEXT)

        when:
        def reader = Files.newInputStream(path)
        then:
        reader.text == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should copy a bucket' () {
        given:
        def folder = Files.createTempDirectory('test')
        def target = folder.resolve('file.txt')
        and:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(50 * 1024)
        final path = s3path("s3://$bucketName/file.txt")
        createObject(path, TEXT)

        when:
        target = FileHelper.copyPath(path, target)

        then:
        target.text == TEXT

        cleanup:
        folder?.deleteDir()
        deleteBucket(bucketName)
    }

    def 'should create a newOutputStream' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(2048)
        final path = s3path("s3://$bucketName/file.txt")

        when:
        def writer = Files.newOutputStream(path)
        TEXT.readLines().each { it ->
            writer.write(it.bytes);
            writer.write((int)('\n' as char))
        }
        writer.close()
        then:
        readObject(path) == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should read a newByteChannel' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(1024)
        final path = s3path("s3://$bucketName/file.txt")
        createObject(path, TEXT)

        when:
        def channel = Files.newByteChannel(path)
        then:
        readChannel(channel, 100) == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should write a byte channel' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(1024)
        final path = s3path("s3://$bucketName/file.txt")

        when:
        def channel = Files.newByteChannel(path, StandardOpenOption.WRITE, StandardOpenOption.CREATE)
        writeChannel(channel, TEXT, 200)
        channel.close()
        then:
        readObject(path) == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should check file size' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(50 * 1024)
        final path = s3path("s3://$bucketName/file.txt")

        when:
        createObject(path, TEXT)
        then:
        Files.size(path) == TEXT.size()

        when:
        Files.size(path.resolve('xxx'))
        then:
        thrown(NoSuchFileException)

        cleanup:
        deleteBucket(bucketName)
    }

    // test fails on GHA, likely due to concurrent execution
    //@IgnoreIf({ System.getenv('GITHUB_ACTIONS') })
    @Ignore
    def 'should list root directory' () {
        given:
        def bucketName1 = createBucket()
        def bucketName2 = createBucket()
        def bucketName3 = createBucket()
        and:
        createObject("$bucketName1/file.1", 'xxx')
        createObject("$bucketName2/foo/file.2", 'xxx')
        createObject("$bucketName2/foo/bar/file.3", 'xxx')

        and:
        def root = s3path('s3://')

        when:
        def paths = Files.newDirectoryStream(root).collect { it.fileName.toString() }
        then:
        paths.contains(bucketName1)
        paths.contains(bucketName2)
        paths.contains(bucketName3)

        when:
        Set<String> dirs = []
        Set<String> files = []
        Files.walkFileTree(root, new SimpleFileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs)
                    throws IOException
            {
                dirs << dir.toString()
                return FileVisitResult.CONTINUE;
            }

            @Override
            FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException
            {
                files << file.toString()
                return FileVisitResult.CONTINUE;
            }

        })
        then:
        println dirs
        println files
        dirs.contains('/')
        dirs.contains("/$bucketName1" as String)
        dirs.contains("/$bucketName2" as String)
        dirs.contains("/$bucketName2/foo" as String)
        dirs.contains("/$bucketName2/foo/bar" as String)
        dirs.contains("/$bucketName3" as String)
        files.contains("/$bucketName1/file.1" as String)
        files.contains("/$bucketName2/foo/file.2" as String)
        files.contains("/$bucketName2/foo/bar/file.3" as String)

        cleanup:
        tryDeleteBucket(bucketName1)
        tryDeleteBucket(bucketName2)
        tryDeleteBucket(bucketName3)
    }


    def 'should stream directory content' () {
        given:
        def bucketName = createBucket()
        createObject("$bucketName/foo/file1.txt",'A')
        createObject("$bucketName/foo/file2.txt",'BB')
        createObject("$bucketName/foo/bar/file3.txt",'CCC')
        createObject("$bucketName/foo/bar/baz/file4.txt",'DDDD')
        createObject("$bucketName/foo/bar/file5.txt",'EEEEE')
        createObject("$bucketName/foo/file6.txt",'FFFFFF')

        when:
        def list = Files.newDirectoryStream(s3path("s3://$bucketName")).collect { it.getFileName().toString() }
        then:
        list.size() == 1
        list == [ 'foo' ]

        when:
        list = Files.newDirectoryStream(s3path("s3://$bucketName/foo")).collect { it.getFileName().toString() }
        then:
        list.size() == 4
        list as Set == [ 'file1.txt', 'file2.txt', 'bar', 'file6.txt' ] as Set

        when:
        list = Files.newDirectoryStream(s3path("s3://$bucketName/foo/bar")).collect { it.getFileName().toString() }
        then:
        list.size() == 3
        list as Set == [ 'file3.txt', 'baz', 'file5.txt' ] as Set

        when:
        list = Files.newDirectoryStream(s3path("s3://$bucketName/foo/bar/baz")).collect { it.getFileName().toString() }
        then:
        list.size() == 1
        list  == [ 'file4.txt' ]

        cleanup:
        deleteBucket(bucketName)
    }


    def 'should check walkTree' () {

        given:
        def bucketName = createBucket()
        createObject("$bucketName/foo/file1.txt",'A')
        createObject("$bucketName/foo/file2.txt",'BB')
        createObject("$bucketName/foo/bar/file3.txt",'CCC')
        createObject("$bucketName/foo/bar/baz/file4.txt",'DDDD')
        createObject("$bucketName/foo/bar/file5.txt",'EEEEE')
        createObject("$bucketName/foo/file6.txt",'FFFFFF')

        when:
        List<String> dirs = []
        Map<String,BasicFileAttributes> files = [:]
        def base = s3path("s3://$bucketName")
        Files.walkFileTree(base, new SimpleFileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException
            {
                dirs << base.relativize(dir).toString()
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException
            {
                files[file.getFileName().toString()] = attrs
                return FileVisitResult.CONTINUE;
            }
        })

        then:
        files.size() == 6
        files ['file1.txt'].size() == 1
        files ['file2.txt'].size() == 2
        files ['file3.txt'].size() == 3
        files ['file4.txt'].size() == 4
        files ['file5.txt'].size() == 5
        files ['file6.txt'].size() == 6
        dirs.size() == 4
        dirs.contains("")
        dirs.contains('foo')
        dirs.contains('foo/bar')
        dirs.contains('foo/bar/baz')


        when:
        dirs = []
        files = [:]
        base = s3path("s3://$bucketName/foo/bar/")
        Files.walkFileTree(base, new SimpleFileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException
            {
                dirs << base.relativize(dir).toString()
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException
            {
                files[file.getFileName().toString()] = attrs
                return FileVisitResult.CONTINUE;
            }
        })

        then:
        files.size()==3
        files.containsKey('file3.txt')
        files.containsKey('file4.txt')
        files.containsKey('file5.txt')
        dirs.size() == 2
        dirs.contains("")
        dirs.contains('baz')

        cleanup:
        deleteBucket(bucketName)
    }

    @Ignore // FIXME 
    def 'should handle dir and files having the same name' () {

        given:
        def bucketName = createBucket()
        createObject("$bucketName/foo",'file-1')
        createObject("$bucketName/foo/bar",'file-2')
        createObject("$bucketName/foo/baz",'file-3')
        and:
        def root = s3path("s3://$bucketName")

        when:
        def file1 = root.resolve('foo')
        then:
        Files.isRegularFile(file1)
        !Files.isDirectory(file1)
        file1.text == 'file-1'

        when:
        def dir1 = root.resolve('foo/')
        then:
        !Files.isRegularFile(dir1)
        Files.isDirectory(dir1)

        when:
        def file2 = root.resolve('foo/bar')
        then:
        Files.isRegularFile(file2)
        !Files.isDirectory(file2)
        file2.text == 'file-2'


        when:
        def parent = file2.parent
        then:
        !Files.isRegularFile(parent)
        Files.isDirectory(parent)

        when:
        Set<String> dirs = []
        Map<String,BasicFileAttributes> files = [:]
        Files.walkFileTree(root, new SimpleFileVisitor<Path>() {

            @Override
            FileVisitResult preVisitDirectory(Path dir, BasicFileAttributes attrs) throws IOException
            {
                dirs << root.relativize(dir).toString()
                return FileVisitResult.CONTINUE;
            }

            @Override
            public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException
            {
                files[root.relativize(file).toString()] = attrs
                return FileVisitResult.CONTINUE;
            }
        })
        then:
        dirs.size() == 2
        dirs.contains('')
        dirs.contains('foo')
        files.size() == 3
        files.containsKey('foo')
        files.containsKey('foo/bar')
        files.containsKey('foo/baz')

        cleanup:
        deleteBucket(bucketName)

    }

    def 'should handle file names with same prefix' () {
        given:
        def bucketName = createBucket()
        and:
        createObject("$bucketName/transcript_index.junctions.fa", 'foo')
        createObject("$bucketName/alpha-beta/file1", 'bar')
        createObject("$bucketName/alpha/file2", 'baz')

        expect:
        Files.exists(s3path("s3://$bucketName/transcript_index.junctions.fa"))
        !Files.exists(s3path("s3://$bucketName/transcript_index.junctions"))
        Files.exists(s3path("s3://$bucketName/alpha-beta/file1"))
        Files.exists(s3path("s3://$bucketName/alpha/file2"))
        Files.exists(s3path("s3://$bucketName/alpha-beta/"))
        Files.exists(s3path("s3://$bucketName/alpha-beta"))
        Files.exists(s3path("s3://$bucketName/alpha/"))
        Files.exists(s3path("s3://$bucketName/alpha"))

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should tag a file' () {
        given:
        def bucketName = createBucket()
        and:
        def path = s3path("s3://$bucketName/alpha.txt")
        def copy = s3path("s3://$bucketName/omega.txt")
        and:
        def client = path.getFileSystem().getClient()

        when:
        path.setTags(FOO: 'Hello world', BAR: 'xyz')
        Files.createFile(path)
        then:
        Files.exists(path)
        and:
        def tags = client .getObjectTags(path.getBucket(), path.getKey())
        tags.find { it.key=='FOO' }.value == 'Hello world'
        tags.find { it.key=='BAR' }.value == 'xyz'

        when:
        copy.setTags(FOO: 'Hola mundo', BAZ: '123')
        Files.copy(path, copy)
        then:
        Files.exists(copy)
        and:
        def copyTags = client .getObjectTags(copy.getBucket(), copy.getKey())
        copyTags.find { it.key=='FOO' }.value == 'Hola mundo'
        copyTags.find { it.key=='BAZ' }.value == '123'
        copyTags.find { it.key=='BAR' } == null

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should download file from encrypted bucket' () {
        given:
        def folder = Files.createTempDirectory('test')
        def target = folder.resolve('test-data.txt')
        and:
        def source = s3path("s3://nf-kms-xyz/test-data.txt")

        when:
        FileHelper.copyPath(source, target)
        then:
        target.exists()
        
        cleanup:
        folder?.deleteDir()
    }

    def 'should upload file to encrypted bucket' () {
        given:
        def KEY = 'arn:aws:kms:eu-west-1:195996028523:key/e97ecf28-951e-4700-bf22-1bd416ec519f'
        and:
        def folder = Files.createTempDirectory('test')
        def source = folder.resolve('hello.txt'); source.text = 'Hello world'
        and:
        def target = s3path("s3://nf-kms-xyz/test-${UUID.randomUUID()}.txt")
        and: // assign some tags
        target.setTags([ONE: 'HELLO'])

        when:
        FileHelper.copyPath(source, target)
        then:
        target.exists()
        
        expect:
        target.getFileSystem().getClient().getObjectKmsKeyId(target.bucket, target.key) == KEY
        and:
        target.getFileSystem().getClient().getObjectTags(target.bucket, target.key).find { it.key=='ONE' }.value == 'HELLO'

        cleanup:
        Files.deleteIfExists(target)
        folder?.deleteDir()
    }

    def 'should upload directory to encrypted bucket' () {
        given:
        def KEY = 'arn:aws:kms:eu-west-1:195996028523:key/e97ecf28-951e-4700-bf22-1bd416ec519f'
        and:
        def folder = Files.createTempDirectory('test')
        def source = folder.resolve('data'); source.mkdir()
        source.resolve('file-1.txt').text = 'file 1'
        source.resolve('file-2.txt').text = 'file 2'
        source.resolve('alpha').mkdir()
        source.resolve('alpha/file-3.txt').text = 'file 3'
        source.resolve('alpha/beta').mkdir()
        source.resolve('alpha/beta/file-4.txt').text = 'file 4'
        source.resolve('alpha/beta/file-5.txt').text = 'file 5'

        and:
        def target = s3path("s3://nf-kms-xyz/test-${UUID.randomUUID()}")
        and: // assign some tags
        target.setTags([ONE: 'HELLO'])

        when:
        FileHelper.copyPath(source, target)
        then:
        target.exists()
        target.resolve('file-1.txt').text == 'file 1'
        target.resolve('file-2.txt').text == 'file 2'
        target.resolve('alpha/file-3.txt').text == 'file 3'
        target.resolve('alpha/beta/file-4.txt').text == 'file 4'
        target.resolve('alpha/beta/file-5.txt').text == 'file 5'

        expect:
        def client = target.getFileSystem().getClient()
        and:
        client.getObjectKmsKeyId(target.bucket,  "$target.key/file-1.txt") == KEY
        client.getObjectKmsKeyId(target.bucket,  "$target.key/alpha/beta/file-5.txt") == KEY
        and:
        client.getObjectTags(target.bucket,  "$target.key/file-1.txt") == [ new Tag('ONE','HELLO') ]
        client.getObjectTags(target.bucket,  "$target.key/alpha/beta/file-5.txt") == [ new Tag('ONE','HELLO') ]

        cleanup:
        target?.deleteDir()
        folder?.deleteDir()
    }

    def 'should download s3 dir to local dir' () {
        given:
        def bucketName = createBucket()
        createObject("$bucketName/cache/foo/file-1",'File one')
        createObject("$bucketName/cache/foo/bar/file-2",'File two')
        createObject("$bucketName/cache/foo/baz/file-3",'File three')
        and:
        def local = Files.createTempDirectory('test')
        def cache1 = local.resolve('cache1')
        def cache2 = local.resolve('cache2')
        and:
        def remote = s3path("s3://$bucketName/cache/foo")

        when:
        CopyMoveHelper.copyToForeignTarget(remote, cache1)
        then:
        cache1.resolve('file-1').exists()
        cache1.resolve('bar/file-2').exists()
        cache1.resolve('baz/file-3').exists()

        when:
        // the use of 'FileHelper.copyPath' will invoke the s3 provider download directory method
        // make sure the resulting local directory structure matches the one created by 'CopyMoveHelper.copyToForeignTarget'
        FileHelper.copyPath(remote, cache2)
        then:
        cache2.resolve('file-1').exists()
        cache2.resolve('bar/file-2').exists()
        cache2.resolve('baz/file-3').exists()

        cleanup:
        local?.deleteDir()
        deleteBucket(bucketName)
    }


    def 'should upload local dir to s3 directory' () {
        given:
        def bucketName = createBucket()
        def local = Files.createTempDirectory('test')
        local.resolve('cache/foo').mkdirs()
        local.resolve('cache/foo/bar').mkdirs()
        local.resolve('cache/foo/baz').mkdirs()
        and:
        local.resolve('cache/foo/file-1').text = 'File one'
        local.resolve('cache/foo/bar/file-2').text = 'File two'
        local.resolve('cache/foo/baz/file-3').text = 'File three'
        local.resolve('cache/foo/baz/file-4').text = 'File four'

        when:
        CopyMoveHelper.copyToForeignTarget(local.resolve('cache/foo'), s3path("s3://$bucketName/cache1"))
        then:
        Files.exists(s3path("s3://$bucketName/cache1/file-1"))
        Files.exists(s3path("s3://$bucketName/cache1/bar/file-2"))
        Files.exists(s3path("s3://$bucketName/cache1/baz/file-3"))
        Files.exists(s3path("s3://$bucketName/cache1/baz/file-4"))

        when:
        FileHelper.copyPath(local.resolve('cache/foo'), s3path("s3://$bucketName/cache2"))
        then:
        Files.exists(s3path("s3://$bucketName/cache2/file-1"))
        Files.exists(s3path("s3://$bucketName/cache2/bar/file-2"))
        Files.exists(s3path("s3://$bucketName/cache2/baz/file-3"))
        Files.exists(s3path("s3://$bucketName/cache2/baz/file-4"))

        cleanup:
        local?.deleteDir()
        deleteBucket(bucketName)
    }

    void "should upload a stream with multiple flush"(){
        given:
        def bucketName = createBucket()
        and:
        def path = s3path("s3://$bucketName/alpha.txt")

        when:
        PrintWriter writer = new PrintWriter(Files.newBufferedWriter(path, Charset.defaultCharset()))
        writer.println '*'*20
        writer.flush()
        writer.println '*'*20
        writer.flush()
        writer.close()

        then:
        Files.readString(s3path("s3://$bucketName/alpha.txt")).length() == 42 // 2*20 + 2 return lines

        cleanup:
        deleteBucket(bucketName)
    }

    void "should upload a stream without flush"(){
        given:
        def bucketName = createBucket()
        and:
        def path = s3path("s3://$bucketName/alpha.txt")

        when:
        PrintWriter writer = new PrintWriter(Files.newBufferedWriter(path, Charset.defaultCharset()))
        writer.println '*'*20
        writer.println '*'*20
        writer.close()

        then:
        Files.readString(s3path("s3://$bucketName/alpha.txt")).length() == 42 // 2*20 + 2 return lines

        cleanup:
        deleteBucket(bucketName)
    }

    @Unroll
    def 'should upload, copy and download a file' () {
        given:
        def TEXT = randomText(FILE_SIZE)
        def folder = Files.createTempDirectory('test')
        def file = Files.write(folder.resolve('foo.data'), TEXT.bytes)
        and:
        def bucket1 = createBucket()
        def bucket2 = createBucket()

        // upload a file to a remote bucket
        when:
        def target1 = s3path("s3://$bucket1/foo.data")
        FileHelper.copyPath(file, target1)
        // the file exist
        then:
        Files.exists(target1)
        Files.size(target1) == Files.size(file)

        // copy a file across buckets
        when:
        def target2 = s3path("s3://$bucket2/foo.data")
        FileHelper.copyPath(target1, target2)
        // the file exist
        then:
        Files.exists(target2)
        Files.size(target2) == Files.size(target1)

        // download a file locally
        when:
        def result = folder.resolve('result.data')
        FileHelper.copyPath(target2, result)
        then:
        Files.exists(result)
        and:
        Files.size(target2) == Files.size(result)

        cleanup:
        deleteBucket(bucket1)
        deleteBucket(bucket2)
        folder?.deleteDir()

        // check the limits in the file `amazon.properties`
        // in the test resources
        where:
        _ | FILE_SIZE
        _ | 50 * 1024
        _ | 11 * 1024 * 1024
    }

    @Unroll
    def 'should set file media type' () {
        given:
        def TEXT = randomText(FILE_SIZE)
        def folder = Files.createTempDirectory('test')
        def file = Files.write(folder.resolve('foo.data'), TEXT.bytes)
        and:
        def bucket1 = createBucket()
        def bucket2 = createBucket()

        // upload a file to a remote bucket
        when:
        def target1 = s3path("s3://$bucket1/foo.data")
        and:
        target1.setContentType('text/foo')
        def client = target1.getFileSystem().getClient()
        and:
        FileHelper.copyPath(file, target1)
        // the file exist
        then:
        Files.exists(target1)
        and:
        client
                .getObjectMetadata(target1.getBucket(), target1.getKey())
                .getContentType() == 'text/foo'

        // copy a file across buckets
        when:
        def target2 = s3path("s3://$bucket2/foo.data")
        and:
        target2.setContentType('text/bar')
        and:
        FileHelper.copyPath(target1, target2)
        // the file exist
        then:
        Files.exists(target2)
        client
                .getObjectMetadata(target2.getBucket(), target2.getKey())
                .getContentType() == 'text/bar'

        cleanup:
        deleteBucket(bucket1)
        deleteBucket(bucket2)
        folder?.deleteDir()

        // check the limits in the file `amazon.properties`
        // in the test resources
        where:
        _ | FILE_SIZE
        _ | 50 * 1024
        _ | 11 * 1024 * 1024
    }

    @Unroll
    def 'should set file storage class' () {
        given:
        def TEXT = randomText(FILE_SIZE)
        def folder = Files.createTempDirectory('test')
        def file = Files.write(folder.resolve('foo.data'), TEXT.bytes)
        and:
        def bucket1 = createBucket()
        def bucket2 = createBucket()

        // upload a file to a remote bucket
        when:
        def target1 = s3path("s3://$bucket1/foo.data")
        and:
        target1.setStorageClass('REDUCED_REDUNDANCY')
        def client = target1.getFileSystem().getClient()
        and:
        FileHelper.copyPath(file, target1)
        // the file exist
        then:
        Files.exists(target1)
        and:
        client
                .getObjectMetadata(target1.getBucket(), target1.getKey())
                .getStorageClass() == 'REDUCED_REDUNDANCY'

        // copy a file across buckets
        when:
        def target2 = s3path("s3://$bucket2/foo.data")
        and:
        target2.setStorageClass('STANDARD_IA')
        and:
        FileHelper.copyPath(target1, target2)
        // the file exist
        then:
        Files.exists(target2)
        client
                .getObjectMetadata(target2.getBucket(), target2.getKey())
                .getStorageClass() == 'STANDARD_IA'

        cleanup:
        deleteBucket(bucket1)
        deleteBucket(bucket2)
        folder?.deleteDir()

        // check the limits in the file `amazon.properties`
        // in the test resources
        where:
        _ | FILE_SIZE
        _ | 50 * 1024
        _ | 11 * 1024 * 1024
    }

    def 'should overwrite a file' () {
        given:
        def bucket1 = createBucket()
        def path = s3path("s3://$bucket1/foo/bar.txt")
        and:
        path.text = 'foo'

        when:
        def file = TraceHelper.newFileWriter(path, true, 'Test')
        file.write('Hola')
        file.close()
        then:
        path.text == 'Hola'

        cleanup:
        deleteBucket(bucket1)
    }

    def 'should not overwrite a file' () {
        given:
        def bucket1 = createBucket()
        def path = s3path("s3://$bucket1/foo/bar.txt")
        and:
        path.text = 'foo'

        when:
        TraceHelper.newFileWriter(path, false, 'Test')
        then:
        def e = thrown(AbortOperationException)
        e.message == "Test file already exists: ${path.toUriString()} -- enable the 'test.overwrite' option in your config file to overwrite existing files"

        cleanup:
        deleteBucket(bucket1)
    }

}
