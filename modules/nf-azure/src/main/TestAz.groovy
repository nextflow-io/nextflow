import java.nio.file.FileSystem
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path

import com.azure.storage.blob.nio.AzureFileSystem

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestAz {

    static void main(String... args) {
        Map<String, Object> config = new HashMap<>();
//        config.put(AzureFileSystem.AZURE_STORAGE_ACCOUNT_KEY, "/p/AB/MIKOjU7y1vmt7lF3f0c3Ubztc8Js+9eKZJyCo58Zpb7y+CLCHOgn/5WR1F2iMV9i12/Yy8g1NAhd74Wg==");
        config.put(AzureFileSystem.AZURE_STORAGE_ACCOUNT_KEY, "h3cUlUTlzgFecGm+/jHrB/AIPdGeuaupYKWIFmjsUWFihOw9RbW5rnKd2/hnd+Jg4ot1YccDtk7yhpm84fkMsQ==");
        config.put(AzureFileSystem.AZURE_STORAGE_FILE_STORES, "my-data");
        FileSystem myFs = FileSystems.newFileSystem(new URI("azb://my-data:?account=nfbucket"), config)

        def root = myFs.getRootDirectories().iterator().next()

        for (Path p : Files.newDirectoryStream(myFs.getPath('my-data:'))) {
            System.out.println(p.toString())
        }


        myFs.getPath('foo').resolve('bar').text = 'Ciao'
    }


}
