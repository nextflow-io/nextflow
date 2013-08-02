package nextflow.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DxHelperTest extends Specification {

    def testDownload() {

        when:
        File target = File.createTempFile('dx-download',null)
        DxHelper.downloadFile( 'file-B7vQg7j0j583xBb7QJV00284', target )

        then:
        target.text.startsWith('>1aboA')

        cleanup:
        target?.delete()

    }
}
