package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires( { System.getenv('AZURE_STORAGE_ACCOUNT_KEY') } )
class AzPathFactoryTest extends Specification {

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

    @Unroll
    def 'should convert az path to uri string' () {
        expect:
        FileHelper.asPath(STR).toUriString() == EXPECTED

        where:
        STR                             | EXPECTED
        'azb://my-data:/file.txt'       | 'azb://my-data:/file.txt'
    }
}
