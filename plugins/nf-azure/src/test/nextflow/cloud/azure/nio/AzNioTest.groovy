package nextflow.cloud.azure.nio

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
import java.nio.file.attribute.BasicFileAttributeView
import java.nio.file.attribute.BasicFileAttributes

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.BlobServiceClientBuilder
import com.azure.storage.common.StorageSharedKeyCredential
import nextflow.exception.AbortOperationException
import nextflow.trace.TraceHelper
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class AzNioTest extends Specification implements AzBaseSpec {

    @Shared
    BlobServiceClient storageClient

    def setupSpec() {
        def accountName = System.getenv('AZURE_STORAGE_ACCOUNT_NAME')
        def accountKey = System.getenv('AZURE_STORAGE_ACCOUNT_KEY')
        def credential = new StorageSharedKeyCredential(accountName, accountKey);
        def  endpoint = String.format(Locale.ROOT, "https://%s.blob.core.windows.net", accountName);
        storageClient = new BlobServiceClientBuilder().endpoint(endpoint).credential(credential).buildClient();
    }


    def 'should create a blob' () {
        given:
        def bucket = createBucket()
        def path = Paths.get(new URI("az://$bucket/file-name.txt"))

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
        def path = Paths.get(new URI("az://$bucket/file-name.txt"))

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
        def path = Paths.get(new URI("az://$bucket/file-name.txt"))
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
        def keyName = "$bucketName/data/alpha.txt"
        createObject(keyName, TEXT)
        and:

        //
        // -- readAttributes
        //
        def path = Paths.get(new URI("az://$keyName"))
        def attrs = Files.readAttributes(path, BasicFileAttributes)
        then:
        attrs.isRegularFile()
        !attrs.isDirectory()
        attrs.size() == 12
        !attrs.isSymbolicLink()
        !attrs.isOther()
        attrs.fileKey() == "/$keyName"
        attrs.lastAccessTime() == null
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
        when:
        def view = Files.getFileAttributeView(path, BasicFileAttributeView)
        then:
        view.readAttributes() == attrs

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
        attrs.fileKey() == "/$bucketName/data/"
        attrs.lastAccessTime() == null
        attrs.lastModifiedTime() == null
        attrs.creationTime() == null

        //
        // -- readAttributes for a bucket
        //
        when:
        attrs = Files.readAttributes(Paths.get(new URI("az://$bucketName")), BasicFileAttributes)
        then:
        !attrs.isRegularFile()
        attrs.isDirectory()
        attrs.size() == 0
        !attrs.isSymbolicLink()
        !attrs.isOther()
        attrs.fileKey() == "/$bucketName"
        attrs.creationTime() == null
        attrs.lastAccessTime() == null
        attrs.lastModifiedTime().toMillis()-start < 5_000

        cleanup:
        if( bucketName ) deleteBucket(bucketName)
    }


    def 'should copy a stream to bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = Paths.get(new URI("az://$bucketName/data/file.txt"))
        
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
        def target = Paths.get(new URI("az://$bucketName/data/file.txt"))
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
        def target = Paths.get(new URI("az://$bucketName/target/file.txt"))

        and:
        final objectName = "$bucketName/source/file.txt"
        final source = Paths.get(new URI("az://$objectName"))
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

    def 'move a remote file to a bucket' () {
        given:
        def TEXT = "Hello world!"

        when:
        def bucketName = createBucket()
        def target = Paths.get(new URI("az://$bucketName/target/file.txt"))

        and:
        final objectName = "$bucketName/source/file.txt"
        final source = Paths.get(new URI("az://$objectName"))
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

    def 'should create a directory' () {

        given:
        def bucketName = getRndBucketName()
        def dir = Paths.get(new URI("az://$bucketName"))

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
        def dir = Paths.get(new URI("az://$bucketName/alpha/bravo/omega/"))
        when:
        Files.createDirectories(dir)
        then:
        Files.exists(Paths.get(new URI("az://$bucketName/alpha/")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha/bravo/")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha/bravo/omega/")))

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
        def path = Paths.get(new URI("az://$bucketName/data/file.txt"))
        Files.createFile(path)
        then:
        existsPath(path)

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create temp file and directory' () {
        given:
        def bucketName = createBucket()
        def base = Paths.get(new URI("az://$bucketName"))

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
        def target = Paths.get(new URI("az://$bucketName/data/file.txt"))
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

    def 'should delete a container' () {
        given:
        final bucketName = createBucket()

        when:
        Files.delete(Paths.get(new URI("az://$bucketName")))
        then:
        !existsPath(bucketName)

    }

    def 'should throw when deleting a not empty container' () {
        given:
        def bucketName = createBucket()
        and:
        createObject("$bucketName/this/that", 'HELLO')

        when:
        def path1 = new URI("az://$bucketName")
        Files.delete(Paths.get(path1))
        then:
        thrown(DirectoryNotEmptyException)

        when:
        def path2 = new URI("az://$bucketName/this")
        Files.delete(Paths.get(path2))
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

    def 'should throw a NoSuchFileException when deleting an object not existing' () {

        given:
        def bucketName = getRndBucketName()
        def path = Paths.get(new URI("az://$bucketName/alpha/bravo"))

        when:
        Files.delete(path)
        then:
        thrown(NoSuchFileException)

    }

    def 'should validate exists method' () {
        given:
        def bucketName = createBucket()
        and:
        def missingBucket = getRndBucketName()
        and:
        createObject("$bucketName/file.txt", 'HELLO')

        expect:
        Files.exists(Paths.get(new URI("az://$bucketName")))
        Files.exists(Paths.get(new URI("az://$bucketName/file.txt")))
        !Files.exists(Paths.get(new URI("az://$bucketName/fooooo.txt")))
        !Files.exists(Paths.get(new URI("az://$missingBucket")))

        cleanup:
        deleteBucket(bucketName)
    }


    def 'should check is it is a directory' () {
        given:
        def bucketName = createBucket()

        when:
        def path = Paths.get(new URI("az://$bucketName"))
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
        def file1 = Paths.get(new URI("az://some/data/file.txt"))
        def file2 = Paths.get(new URI("az://some/data/file.txt"))
        def file3 = Paths.get(new URI("az://some/data/fooo.txt"))

        expect:
        Files.isSameFile(file1, file2)
        !Files.isSameFile(file1, file3)

    }

    def 'should create a newBufferedReader' () {
        given:
        def bucketName = createBucket()
        and:
        def TEXT = randomText(50 * 1024)
        def path = Paths.get(new URI("az://$bucketName/file.txt"))
        createObject(path, TEXT)

        when:
        def reader = Files.newBufferedReader(path, Charset.forName('UTF-8'))
        then:
        reader.text == TEXT

        when:
        def unknown = Paths.get(new URI("az://$bucketName/unknown.txt"))
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
        final path = Paths.get(new URI("az://$bucketName/file.txt"))

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
        final path = Paths.get(new URI("az://$bucketName/file.txt"))
        createObject(path, TEXT)

        when:
        def reader = Files.newInputStream(path)
        then:
        reader.text == TEXT

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should create a newOutputStream' () {
        given:
        def bucketName = createBucket()
        and:
        final TEXT = randomText(2048)
        final path = Paths.get(new URI("az://$bucketName/file.txt"))

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
        final path = Paths.get(new URI("az://$bucketName/file.txt"))
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
        final path = Paths.get(new URI("az://$bucketName/file.txt"))

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
        final path = Paths.get(new URI("az://$bucketName/file.txt"))

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
    @IgnoreIf({ System.getenv('GITHUB_ACTIONS') })
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
        def root = Paths.get(new URI('az:///'))

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
        def list = Files.newDirectoryStream(Paths.get(new URI("az://$bucketName"))).collect { it.getFileName().toString() }
        then:
        list.size() == 1
        list == [ 'foo' ]

        when:
        list = Files.newDirectoryStream(Paths.get(new URI("az://$bucketName/foo"))).collect { it.getFileName().toString() }
        then:
        list.size() == 4
        list as Set == [ 'file1.txt', 'file2.txt', 'bar', 'file6.txt' ] as Set

        when:
        list = Files.newDirectoryStream(Paths.get(new URI("az://$bucketName/foo/bar"))).collect { it.getFileName().toString() }
        then:
        list.size() == 3
        list as Set == [ 'file3.txt', 'baz', 'file5.txt' ] as Set

        when:
        list = Files.newDirectoryStream(Paths.get(new URI("az://$bucketName/foo/bar/baz"))).collect { it.getFileName().toString() }
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
        def base = Paths.get(new URI("az://$bucketName"))
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
        base = Paths.get(new URI("az://$bucketName/foo/bar/"))
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

    def 'should handle dir and files having the same name' () {

        given:
        def bucketName = createBucket()
        createObject("$bucketName/foo",'file-1')
        createObject("$bucketName/foo/bar",'file-2')
        createObject("$bucketName/foo/baz",'file-3')
        and:
        def root = Paths.get(new URI("az://$bucketName"))

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
        Files.exists(Paths.get(new URI("az://$bucketName/transcript_index.junctions.fa")))
        !Files.exists(Paths.get(new URI("az://$bucketName/transcript_index.junctions")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha-beta/file1")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha/file2")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha-beta/")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha-beta")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha/")))
        Files.exists(Paths.get(new URI("az://$bucketName/alpha")))

        cleanup:
        deleteBucket(bucketName)
    }

    def 'should overwrite a file' () {
        given:
        def bucket1 = createBucket()
        def path = Path.of(new URI("az://$bucket1/foo/bar.txt"))
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
        def path = Path.of(new URI("az://$bucket1/foo/bar.txt"))
        and:
        path.text = 'foo'

        when:
        TraceHelper.newFileWriter(path, false, 'Test')
        then:
        def e = thrown(AbortOperationException)
        e.message == "Test file already exists: ${path.toUriString()}"

        cleanup:
        deleteBucket(bucket1)
    }

}
