package nextflow.util

import groovy.json.JsonSlurper
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


    def testUpload() {

        setup:
        def home = System.getProperty('user.home')
        def config = new File(home, '.dnanexus_config/environment.json')
        def json = new JsonSlurper().parseText(config.text)
        def projectId = json.DX_PROJECT_CONTEXT_ID as String

        when:
        File target = File.createTempFile('dxupload', 'test')
        target.text = "Hello there!\n"
        String id = DxHelper.uploadFile(target, 'test-upload', projectId)

        then:
        id.startsWith('file-')

        cleanup:
        target?.delete()
    }


}
