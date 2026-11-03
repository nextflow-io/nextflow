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

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.script.dsl.Description

/**
 * Model Azure Batch registry config settings from nextflow config file
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
@CompileStatic
class AzRegistryOpts implements ConfigScope {

    @ConfigOption
    @Description("""
        The container registry from which to pull the Docker images (default: `docker.io`).
    """)
    final String server

    @ConfigOption
    @Description("""
        The username to connect to a private container registry.
    """)
    final String userName

    @ConfigOption
    @Description("""
        The password to connect to a private container registry.
    """)
    final String password

    AzRegistryOpts() {
        this(Collections.emptyMap())
    }

    AzRegistryOpts(Map config, Map<String,String> env=SysEnv.get()) {
        assert config!=null
        this.server = config.server ?: 'docker.io'
        this.userName = config.userName ?: env.get('AZURE_REGISTRY_USER_NAME')
        this.password = config.password ?: env.get('AZURE_REGISTRY_PASSWORD')
    }

    boolean isConfigured() {
        if( userName && password )
            return true
        if( !userName && !password )
            return false
        throw new IllegalArgumentException("Invalid Container Registry configuration - Make sure userName and password are set for Container Registry")
    }

}
