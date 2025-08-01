/*
 * Copyright 2022, Seqera Labs
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
import nextflow.SysEnv
import nextflow.cloud.azure.nio.AzFileSystemProvider
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.script.dsl.Description

/**
 * Model Azure Entra (formerly Active Directory) config options
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@CompileStatic
class AzActiveDirectoryOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        The service principal client ID. Defaults to environment variable `AZURE_CLIENT_ID`.
    """)
    final String servicePrincipalId

    @ConfigOption
    @Description("""
        The service principal client secret. Defaults to environment variable `AZURE_CLIENT_SECRET`.
    """)
    final String servicePrincipalSecret

    @ConfigOption
    @Description("""
        The Azure tenant ID. Defaults to environment variable `AZURE_TENANT_ID`.
    """)
    final String tenantId

    AzActiveDirectoryOpts(Map config, Map<String, String> env = null) {
        assert config != null
        this.servicePrincipalId = config.servicePrincipalId ?: SysEnv.get('AZURE_CLIENT_ID')
        this.servicePrincipalSecret = config.servicePrincipalSecret ?: SysEnv.get('AZURE_CLIENT_SECRET')
        this.tenantId = config.tenantId ?: SysEnv.get('AZURE_TENANT_ID')
    }

    Map<String, Object> getEnv() {
        Map<String, Object> props = new HashMap<>();
        props.put(AzFileSystemProvider.AZURE_CLIENT_ID, servicePrincipalId)
        props.put(AzFileSystemProvider.AZURE_CLIENT_SECRET, servicePrincipalSecret)
        props.put(AzFileSystemProvider.AZURE_TENANT_ID, tenantId)
        return props
    }

    boolean isConfigured() {
        if (servicePrincipalId && servicePrincipalSecret && tenantId)
            return true
        if (!servicePrincipalId && !servicePrincipalSecret && !tenantId)
            return false
        throw new IllegalArgumentException("Invalid Service Principal configuration - Make sure servicePrincipalId and servicePrincipalClient are set in nextflow.config or configured via environment variables")
    }

}
