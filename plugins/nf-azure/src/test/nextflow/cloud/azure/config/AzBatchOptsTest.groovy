package nextflow.cloud.azure.config

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchOptsTest extends Specification {

    @Unroll
    def 'should fetch location' ()  {

        expect:
        new AzBatchOpts(endpoint: ENDPOINT, location: LOC).getLocation() == EXPECTED

        where:
        LOC                 | ENDPOINT                                  |  EXPECTED
        null                | null                                      | null
        'europe1'           | null                                      | 'europe1'
        'europe1'           | 'https://foo.xyz.batch.azure.com'         | 'europe1'
        null                | 'https://foo.india2.batch.azure.com'      | 'india2'

    }

    @Unroll
    def 'should fetch account name' ()  {

        expect:
        new AzBatchOpts(endpoint: ENDPOINT, accountName: NAME).getAccountName() == EXPECTED

        where:
        NAME                | ENDPOINT |  EXPECTED
        null                | null                                      | null
        'foo'               | null                                      | 'foo'
        'foo'               | 'https://bar.eu1.batch.azure.com'         | 'foo'
        null                | 'https://bar.eu2.batch.azure.com'         | 'bar'

    }
}
