package nextflow.cloud.azure.config

import spock.lang.Specification

/**
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class AzRegistryOptsTest extends Specification {

    def 'should get server, user name & password'() {
        when:
        def opts1 = new AzRegistryOpts([:], [:])
        then:
        opts1.server == 'docker.io'
        opts1.userName == null
        opts1.password == null

        when:
        def opts2 = new AzRegistryOpts(
                [server: 'xyz.io', userName: 'xyz', password: 'A1B2C3'],
                [AZURE_REGISTRY_USER_NAME: 'env-userName', AZURE_REGISTRY_PASSWORD:'env-password'])
        then:
        opts2.server == 'xyz.io'
        opts2.userName == 'xyz'
        opts2.password == 'A1B2C3'


        when:
        def opts3 = new AzRegistryOpts(
                [:],
                [AZURE_REGISTRY_USER_NAME: 'env-userName', AZURE_REGISTRY_PASSWORD:'env-password'])
        then:
        opts3.server == 'docker.io'
        opts3.userName == 'env-userName'
        opts3.password == 'env-password'
    }

}
