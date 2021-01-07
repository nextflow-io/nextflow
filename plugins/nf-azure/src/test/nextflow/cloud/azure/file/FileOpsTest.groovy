package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires( { System.getenv('AZURE_STORAGE_ACCOUNT_KEY') } )
class FileOpsTest extends Specification {

    def setup() {

        def KEY = System.getenv('AZURE_STORAGE_ACCOUNT_KEY')
        def ACCOUNT = "nfbucket"
        def STORES = "my-data"

        Global.session = Mock(Session) {
            getConfig() >> [azure: [storage: [
                    accountKey: KEY,
                    accountName: ACCOUNT,
                    fileStores: STORES
            ]]]
        }
    }

    def 'should create subdirs' () {
        given:
        def directory = FileHelper.asPath('azb://my-data:/work')

        when:
        directory.mkdir()
        then:
        directory.exists()
        directory.isDirectory()

    }

}
