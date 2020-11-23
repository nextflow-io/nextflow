package nextflow.cloud.azure.util

import java.nio.file.Path

import com.azure.storage.blob.nio.AzureFileSystem
import groovy.transform.CompileStatic
import nextflow.file.FileHelper
import nextflow.file.FileSystemPathFactory

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzurePathFactory extends FileSystemPathFactory {

    @Override
    protected Path parseUri(String uri) {
        if( !uri.startsWith('azb://') )
            return null

        Map<String, Object> config = new HashMap<>();
        config.put(AzureFileSystem.AZURE_STORAGE_ACCOUNT_KEY, "/p/AB/MIKOjU7y1vmt7lF3f0c3Ubztc8Js+9eKZJyCo58Zpb7y+CLCHOgn/5WR1F2iMV9i12/Yy8g1NAhd74Wg==");
        config.put(AzureFileSystem.AZURE_STORAGE_FILE_STORES, "my-data");

        final fs = FileHelper.getOrCreateFileSystemFor(new URI(uri), config)
        final str = uri.substring(6)
        final p = str.indexOf(':')
        if( p==-1 )
            throw new IllegalArgumentException("Missing Azure blob container name uri: $uri")
        final container = str.substring(0,p+1)
        final path = stripParams(str.substring(p+1))

        return path ? fs.getPath(container, path) : fs.getPath(container)
    }

    String stripParams(String str) {
        final p = str.indexOf('?')
        return p>0 ? str.substring(0,p) : null
    }

    @Override
    protected String toUriString(Path path) {
        return null
    }

}
