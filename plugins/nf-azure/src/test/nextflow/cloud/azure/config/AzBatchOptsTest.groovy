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
        new AzBatchOpts(endpoint: ENDPOINT, location: LOC, [:]).getLocation() == EXPECTED

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
        new AzBatchOpts(endpoint: ENDPOINT, accountName: NAME, [:]).getAccountName() == EXPECTED

        where:
        NAME                | ENDPOINT                                  |  EXPECTED
        null                | null                                      | null
        'foo'               | null                                      | 'foo'
        'foo'               | 'https://bar.eu1.batch.azure.com'         | 'foo'
        null                | 'https://bar.eu2.batch.azure.com'         | 'bar'

    }

    def 'should get account name & key'() {
        when:
        def opts1 = new AzBatchOpts([:], [:])
        then:
        opts1.accountName == null
        opts1.accountKey == null

        when:
        def opts2 = new AzBatchOpts(
                [accountName: 'xyz', accountKey: '123'],
                [AZURE_BATCH_ACCOUNT_NAME: 'env-name', AZURE_BATCH_ACCOUNT_KEY:'env-key'])
        then:
        opts2.accountName == 'xyz'
        opts2.accountKey == '123'


        when:
        def opts3 = new AzBatchOpts(
                [:],
                [AZURE_BATCH_ACCOUNT_NAME: 'env-name', AZURE_BATCH_ACCOUNT_KEY:'env-key'])
        then:
        opts3.accountName == 'env-name'
        opts3.accountKey == 'env-key'
    }

    @Unroll
    def 'should check azcopy install' () {
        given:
        AzBatchOpts opts = Spy(AzBatchOpts, constructorArgs: [ CONFIG,[:] ])

        expect:
        opts.getCopyToolInstallMode() == EXPECTED

        where:
        EXPECTED                    | CONFIG
        CopyToolInstallMode.task    | [:]
        CopyToolInstallMode.node    | [allowPoolCreation: true]
        CopyToolInstallMode.node    | [autoPoolMode: true]
        CopyToolInstallMode.node    | [allowPoolCreation: true, copyToolInstallMode: 'node']
        CopyToolInstallMode.task    | [allowPoolCreation: true, copyToolInstallMode: 'task']
        CopyToolInstallMode.task    | [copyToolInstallMode: 'task']
        CopyToolInstallMode.node    | [copyToolInstallMode: 'node']

    }
}
