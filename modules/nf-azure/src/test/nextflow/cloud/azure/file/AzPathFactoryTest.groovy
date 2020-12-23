package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.file.FileHelper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzPathFactoryTest extends Specification {

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

    @Unroll
    def 'should convert az path to uri string' () {
        expect:
        FileHelper.asPath(STR).toUriString() == EXPECTED

        where:
        STR                             | EXPECTED
        'azb://my-data:/file.txt'       | 'azb://my-data:/file.txt'
    }
}
