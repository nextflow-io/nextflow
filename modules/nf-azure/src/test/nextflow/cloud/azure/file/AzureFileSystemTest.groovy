package nextflow.cloud.azure.file

import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import com.azure.storage.blob.nio.AzureFileSystem
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzureFileSystemTest extends Specification {

    @Shared
    AzureFileSystem fileSystem

    String MY_KEY = "h3cUlUTlzgFecGm+/jHrB/AIPdGeuaupYKWIFmjsUWFihOw9RbW5rnKd2/hnd+Jg4ot1YccDtk7yhpm84fkMsQ=="
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

    def 'should create directory' () {
        given:
        def dir = fileSystem.getPath( 'mydir')
        when:
        Files.createDirectory(dir)
        then:
        Files.exists(dir)               
        Files.isDirectory(dir)
    }

}
