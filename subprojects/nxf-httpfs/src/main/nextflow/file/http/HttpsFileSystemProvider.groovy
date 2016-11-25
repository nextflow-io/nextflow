package nextflow.file.http

import groovy.transform.CompileStatic

/**
 * JSR-203 file system provider for HTTPS protocol
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class HttpsFileSystemProvider extends XFileSystemProvider {

    @Override
    String getScheme() {
        return 'https'
    }

}
