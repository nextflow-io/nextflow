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

        def KEY = "M47dC/OYTG40e76tuUSqPtcAAbhAmnL1g8auiTHOVra0E3JKw+EfsPc3jegqPNsPene9CmGblvnai/SEU8Oa2w=="
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
