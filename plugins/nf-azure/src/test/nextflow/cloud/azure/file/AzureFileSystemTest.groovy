package nextflow.cloud.azure.file

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import com.azure.storage.blob.nio.AzureFileSystem
import spock.lang.Requires
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires( { System.getenv('AZURE_STORAGE_ACCOUNT_KEY') } )
class AzureFileSystemTest extends Specification {

    @Shared
    AzureFileSystem fileSystem

    String MY_KEY = System.getenv('AZURE_STORAGE_ACCOUNT_KEY')
    String MY_ACCOUNT = "nfbucket"
    String MY_CONTAINER = "my-data"

    def setup() {
        Map<String, Object> config = new HashMap<>();
        config.put(AzureFileSystem.AZURE_STORAGE_ACCOUNT_KEY, MY_KEY);
        config.put(AzureFileSystem.AZURE_STORAGE_FILE_STORES, MY_CONTAINER);
        fileSystem = (AzureFileSystem) FileSystems.newFileSystem(new URI("azb://?account=$MY_ACCOUNT"), config);
    }
            

    def 'should list files' () {

        when:
        def dirPath = fileSystem.getPath("${MY_CONTAINER}:")
        def list = []
        for (Path p : Files.newDirectoryStream(dirPath)) {
            list.add( p.toString() )
        }

        println list
        then:
        list.size()
    }

    def 'should get a file' () {
        when:
        def path1 = fileSystem.getPath("${MY_CONTAINER}:", '/file1.txt')
        def text1 = Files.readAllLines(path1)
        println(text1)
        then:
        true

        when:
        def path2 = fileSystem.getPath("${MY_CONTAINER}:/file1.txt")
        def text2 = Files.readAllLines(path2)
        then:
        text1 == text2

    }

    def 'should create empty dir' () {
        when:
        def dir = fileSystem.getPath( 'mydir')
        Files.createDirectory(dir)
        then:
        Files.exists(dir)
        Files.isDirectory(dir)
    }

    def 'should create directory and a file' () {
        given:
        def dir = fileSystem.getPath( 'mydir')
        when:
        Files.createDirectory(dir)
        then:
        Files.exists(dir)               
        Files.isDirectory(dir)

        when:
        def file = dir.resolve('my-file.txt')
        Files.write(file, 'Hello world'.bytes)
        then:
        Files.exists(file)
        Files.readAllLines(file)[0] == 'Hello world'
    }

    def 'should create directories' () {
        when:
        def dirs = fileSystem.getPath( 'mydir/mydir2/mydir3')
        Files.createDirectories(dirs)
        then:
        Files.isDirectory(fileSystem.getPath('mydir'))
        Files.isDirectory(fileSystem.getPath('mydir/mydir2'))
        Files.isDirectory(fileSystem.getPath('mydir/mydir2/mydir3'))
    }

    def 'should create file' () {
        when:
        def file = fileSystem.getPath( 'my-file.txt')
        Files.createFile(file)
        then:
        Files.exists(file)

        /*
         java.lang.UnsupportedOperationException
            at com.azure.storage.blob.nio.AzureFileSystemProvider.newByteChannel(AzureFileSystemProvider.java:271)
            at java.nio.file.Files.newByteChannel(Files.java:361)
            at java.nio.file.Files.createFile(Files.java:632)
         */
    }
}
