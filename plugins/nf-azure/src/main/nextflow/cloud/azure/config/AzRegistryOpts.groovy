/*
 * Copyright 2021, Microsoft Corp
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
 * Model Azure Batch registry config settings from nextflow config file
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@CompileStatic
class AzRegistryOpts {

    private Map<String,String> sysEnv

    String server
    String userName
    String password

    AzRegistryOpts() {
        this(Collections.emptyMap())
    }

    AzRegistryOpts(Map config, Map<String,String> env=null) {
        assert config!=null
        this.sysEnv = env==null ? new HashMap<String,String>(System.getenv()) : env
        this.server = config.server ?: 'docker.io'
        this.userName = config.userName ?: sysEnv.get('AZURE_REGISTRY_USER_NAME')
        this.password = config.password ?: sysEnv.get('AZURE_REGISTRY_PASSWORD')
    }

    boolean isConfigured() {
        return this.userName && this.password
    }

    boolean isIncomplete() {
        //check whether one of the two mandatory options is missing, but not both (which would me, no config is provided)
        return (!this.userName && this.password) || (this.userName && !this.password)
    }
}
