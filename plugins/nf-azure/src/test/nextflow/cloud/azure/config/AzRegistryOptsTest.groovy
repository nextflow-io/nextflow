/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

    def 'should get managed identity resource id'() {
        expect: 'config value is used'
        new AzRegistryOpts([managedIdentityResourceId: 'res-id'], [:]).managedIdentityResourceId == 'res-id'

        and: 'env fallback is used when config absent'
        new AzRegistryOpts([:], [AZURE_REGISTRY_MANAGED_RESOURCE_ID: 'env-res-id']).managedIdentityResourceId == 'env-res-id'

        and: 'config value wins over env'
        new AzRegistryOpts([managedIdentityResourceId: 'res-id'], [AZURE_REGISTRY_MANAGED_RESOURCE_ID: 'env-res-id']).managedIdentityResourceId == 'res-id'

        and: 'null when neither set'
        new AzRegistryOpts([:], [:]).managedIdentityResourceId == null
    }

    def 'should validate isConfigured'() {
        expect:
        new AzRegistryOpts([:], [:]).isConfigured() == false
        new AzRegistryOpts([userName: 'foo', password: 'bar'], [:]).isConfigured() == true
        new AzRegistryOpts([managedIdentityResourceId: 'res-id'], [:]).isConfigured() == true
        new AzRegistryOpts([managedIdentityResourceId: 'res-id'], [:]).usesManagedIdentity() == true
        new AzRegistryOpts([userName: 'foo', password: 'bar'], [:]).usesManagedIdentity() == false

        when: 'managed identity takes precedence over incomplete user/password'
        def opts = new AzRegistryOpts([managedIdentityResourceId: 'res-id', userName: 'foo'], [:])
        then:
        opts.isConfigured() == true
        opts.usesManagedIdentity() == true

        when: 'partial user/password throws'
        new AzRegistryOpts([userName: 'foo'], [:]).isConfigured()
        then:
        thrown(IllegalArgumentException)
    }

}
