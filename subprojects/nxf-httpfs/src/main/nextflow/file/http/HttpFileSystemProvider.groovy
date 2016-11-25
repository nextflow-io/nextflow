package nextflow.file.http

import groovy.transform.CompileStatic

/**
 * JSR-203 file system provider for HTTP protocol
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class HttpFileSystemProvider extends XFileSystemProvider {

    private static final String SCHEME = 'http'

    @Override
    String getScheme() {
        return SCHEME
    }
}
