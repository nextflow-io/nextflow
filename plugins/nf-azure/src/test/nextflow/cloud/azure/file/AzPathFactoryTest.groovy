package nextflow.cloud.azure.file

import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.nio.AzPath
import spock.lang.Requires
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Requires({System.getenv('AZURE_STORAGE_ACCOUNT_NAME') && System.getenv('AZURE_STORAGE_ACCOUNT_KEY')})
class AzPathFactoryTest extends Specification {

    def 'should create az azure path' () {
        given:
        def CONFIG = [azure: [
                storage: [
                        accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                        accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
                ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
        and:

        when:
        def path = AzPathFactory.parse(AZ_URI)
        then:
        path instanceof AzPath
        (path as AzPath).containerName == CONTAINER
        (path as AzPath).blobName() == BLOB

        when:
        def ret = AzPathFactory.getUriString(path)
        then:
        ret == AZ_URI

        cleanup:
        Global.session = null

        where:
        AZ_URI                          | CONTAINER     | BLOB
        'az://my-data/foo/bar'          | 'my-data'     | 'foo/bar'
        'az://my-data/data/*{1,2}.fq.gz'| 'my-data'     | 'data/*{1,2}.fq.gz'
    }


    def 'should create az azure path and remove storage account from the beginning of the path, if present' () {
        given:
        def CONFIG = [azure: [
            storage: [
                accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
            ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
        and:

        when:
        def path = AzPathFactory.parse(AZ_URI)
        then:
        path instanceof AzPath
        (path as AzPath).containerName == CONTAINER
        (path as AzPath).blobName() == BLOB

        cleanup:
        Global.session = null

        where:
        storageAccount                              |  AZ_URI                                                       | CONTAINER             | BLOB
        System.getenv('AZURE_STORAGE_ACCOUNT_NAME') | "az://${storageAccount}.my-data/foo/bar"                      | 'my-data'             | 'foo/bar'
        System.getenv('AZURE_STORAGE_ACCOUNT_NAME') | "az://${storageAccount}.my-data/data/*{1,2}.fq.gz"            | 'my-data'             | 'data/*{1,2}.fq.gz'
        System.getenv('AZURE_STORAGE_ACCOUNT_NAME') | "az://${storageAccount}.my-data/${storageAccount}.bar"        | 'my-data'             | "${storageAccount}.bar"
        System.getenv('AZURE_STORAGE_ACCOUNT_NAME') | "az://${storageAccount}.${storageAccount}./foo/bar"           | "${storageAccount}."  | 'foo/bar'
        System.getenv('AZURE_STORAGE_ACCOUNT_NAME') | "az://my-data/${storageAccount}.txt"                          | 'my-data'             | "${storageAccount}.txt"

    }

    def 'should throw illegal path' () {
        given:
        def CONFIG = [azure: [
                storage: [
                        accountKey: System.getenv('AZURE_STORAGE_ACCOUNT_KEY'),
                        accountName: System.getenv('AZURE_STORAGE_ACCOUNT_NAME'),
                ]
        ]]
        Global.session = Mock(Session) { getConfig() >> CONFIG }
        and:

        when:
        AzPathFactory.parse('az:///fooo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Invalid Azure path URI - make sure the schema prefix does not container more than two slash characters - offending value: az:///fooo'
    }
}
