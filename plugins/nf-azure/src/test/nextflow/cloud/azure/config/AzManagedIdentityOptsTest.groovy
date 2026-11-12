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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzManagedIdentityOptsTest extends Specification {

    def 'should create manage identity opts' () {
        when:
        def opts = new AzManagedIdentityOpts([:])
        then:
        opts.clientId == null
        !opts.system

        when:
        opts = new AzManagedIdentityOpts(clientId: 'foo', system: false)
        then:
        opts.clientId == 'foo'
        !opts.system

        when:
        opts = new AzManagedIdentityOpts(system: true)
        then:
        opts.clientId == null
        opts.system

        when:
        opts = new AzManagedIdentityOpts(system: 'false')
        then:
        opts.clientId == null
        !opts.system

        when:
        opts = new AzManagedIdentityOpts(system: 'true')
        then:
        opts.clientId == null
        opts.system
    }

    def 'should fall back to env variables' () {
        when:
        def opts = new AzManagedIdentityOpts([:], [AZURE_MANAGED_IDENTITY_USER: 'from-env'])
        then:
        opts.clientId == 'from-env'
        !opts.system

        when:
        opts = new AzManagedIdentityOpts([:], [AZURE_MANAGED_IDENTITY_SYSTEM: 'true'])
        then:
        opts.clientId == null
        opts.system
    }

    def 'should give config precedence over env variables' () {
        when: 'config clientId wins over env'
        def opts = new AzManagedIdentityOpts([clientId: 'from-config'], [AZURE_MANAGED_IDENTITY_USER: 'from-env'])
        then:
        opts.clientId == 'from-config'

        when: 'config system=false wins over env system=true'
        opts = new AzManagedIdentityOpts([system: false], [AZURE_MANAGED_IDENTITY_SYSTEM: 'true'])
        then:
        !opts.system
    }

    @Unroll
    def 'should get env' () {
        when:
        def opts = new AzManagedIdentityOpts(OPTS)
        then:
        opts.getEnv() == ENV

        where:
        OPTS                | ENV
        [:]                 | [AZURE_MANAGED_IDENTITY_USER: null, AZURE_MANAGED_IDENTITY_SYSTEM: false]
        [clientId:'foo']    | [AZURE_MANAGED_IDENTITY_USER: 'foo', AZURE_MANAGED_IDENTITY_SYSTEM: false]
        [system:'true']     | [AZURE_MANAGED_IDENTITY_USER: null, AZURE_MANAGED_IDENTITY_SYSTEM: true]
    }

    def 'should check is valid' (){
        expect:
        !new AzManagedIdentityOpts([:]).isConfigured()
        and:
        new AzManagedIdentityOpts([clientId: 'foo']).isConfigured()
        new AzManagedIdentityOpts([system: true]).isConfigured()

        when:
        new AzManagedIdentityOpts([clientId: 'foo',system: true]).isConfigured()
        then:
        thrown(IllegalArgumentException)
    }

}
