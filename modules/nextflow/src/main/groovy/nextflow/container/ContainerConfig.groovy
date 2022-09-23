/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.container

import groovy.transform.CompileStatic

/**
 * Models container engine configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ContainerConfig extends LinkedHashMap {

    private Map<String,String> sysEnv = System.getenv()

    /* required by Kryo deserialization -- do not remove */
    private ContainerConfig() { }

    ContainerConfig(Map config) {
        super(config)
    }

    ContainerConfig(Map config, Map<String,String> env) {
        super(config)
        this.sysEnv = env
    }

    boolean isEnabled() {
        get('enabled')?.toString() == 'true'
    }

    String getEngine() {
        get('engine')
    }

    List<String> getEnvWhitelist() {
        def result = get('envWhitelist')
        if( !result )
            return Collections.emptyList()

        if( result instanceof CharSequence )
            return result.tokenize(',').collect { it.trim() }

        if( result instanceof List )
            return result

        throw new IllegalArgumentException("Not a valid `envWhitelist` argument")
    }

    boolean entrypointOverride() {
        def result = get('entrypointOverride')
        if( result == null )
            result = sysEnv.get('NXF_CONTAINER_ENTRYPOINT_OVERRIDE')
        if( result != null )
            return Boolean.parseBoolean(result.toString())
        return false
    }
}
