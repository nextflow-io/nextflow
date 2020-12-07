package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FileOpsTest extends Specification {

    def setup() {

        def KEY = "h3cUlUTlzgFecGm+/jHrB/AIPdGeuaupYKWIFmjsUWFihOw9RbW5rnKd2/hnd+Jg4ot1YccDtk7yhpm84fkMsQ=="
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
