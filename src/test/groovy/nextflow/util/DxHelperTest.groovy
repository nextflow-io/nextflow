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


    def testUpload() {

        setup:
        String projectId = System.getenv('DX_PROJECT_CONTEXT_ID') ?: 'project-B7fQ9vj0FqXB2z80y5FQ0JGG'

        when:
        File target = new File("fileToUpload")
        target.text = "Esto es una prueba"
        String id = DxHelper.uploadFile(target,target.name)

        then:
        id.startsWith('file-')

        cleanup:
        target?.delete()
    }
}
