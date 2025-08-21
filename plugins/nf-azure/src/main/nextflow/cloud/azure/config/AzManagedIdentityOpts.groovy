/*
 * Copyright 2013-2024, Seqera Labs
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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.cloud.azure.nio.AzFileSystemProvider
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.script.dsl.Description

/**
 * Model Azure managed identity config options
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@ToString(includePackage = false, includeNames = true)
@EqualsAndHashCode
@CompileStatic
class AzManagedIdentityOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        The client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview). Defaults to environment variable `AZURE_MANAGED_IDENTITY_USER`.
    """)
    final String clientId

    @ConfigOption
    @Description("""
        When `true`, use the system-assigned [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) to authenticate Azure resources. Defaults to environment variable `AZURE_MANAGED_IDENTITY_SYSTEM`.
    """)
    final boolean system

    AzManagedIdentityOpts(Map config) {
        assert config != null
        this.clientId = config.clientId
        this.system = Boolean.parseBoolean(config.system as String)
    }

    Map<String, Object> getEnv() {
        Map<String, Object> props = new HashMap<>();
        props.put(AzFileSystemProvider.AZURE_MANAGED_IDENTITY_USER, clientId)
        props.put(AzFileSystemProvider.AZURE_MANAGED_IDENTITY_SYSTEM, system)
        return props
    }

    boolean isConfigured() {
        if( clientId && !system )
            return true
        if( !clientId && system )
            return true
        if( !clientId && !system )
            return false
        throw new IllegalArgumentException("Invalid Managed Identity configuration - Make sure the `clientId` or `system` is set in the nextflow config")
    }

}
