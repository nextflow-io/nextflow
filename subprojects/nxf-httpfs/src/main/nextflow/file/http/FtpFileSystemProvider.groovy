package nextflow.file.http

import groovy.transform.CompileStatic

/**
 * JSR-203 file system provider for FTP protocol
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FtpFileSystemProvider extends XFileSystemProvider {

    @Override
    String getScheme() {
        return 'ftp'
    }
}
