/*
 * Copyright 2021, Microsoft Corp
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

/**
 * Model Azure identity options from nextflow config file
 *
 * @author Abhinav Sharma <abhi18av@outlook.com>
 */
@CompileStatic
class AzIdentityOpts {

    private Map<String, String> sysEnv
    String servicePrincipalId
    String servicePrincipalSecret
    String tenantId

    AzIdentityOpts(Map config, Map<String, String> env = null) {
        assert config != null
        this.sysEnv = env == null ? new HashMap<String, String>(System.getenv()) : env
        this.servicePrincipalId = config.servicePrincipalId ?: sysEnv.get('AZURE_SP_CLIENT_ID')
        this.servicePrincipalSecret = config.servicePrincipalSecret ?: sysEnv.get('AZURE_SP_CLIENT_SECRET')
        this.tenantId = config.tenantId ?: sysEnv.get('AZURE_SP_TENANT_ID')
    }


    boolean isConfigured() {
        if (servicePrincipalId && servicePrincipalSecret && tenantId)
            return true
        if (!servicePrincipalId && !servicePrincipalSecret && !tenantId)
            return false
        throw new IllegalArgumentException("Invalid Service Principal configuration - Make sure servicePrincipalId and servicePrincipalClient are set in nextflow.config or configured via environment variables")
    }

}
