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

        where:
        AZ_URI                          | CONTAINER     | BLOB
        'az://my-data/foo/bar'          | 'my-data'     | 'foo/bar'
        'az://my-data/data/*{1,2}.fq.gz'| 'my-data'     | 'data/*{1,2}.fq.gz'
    }

}
