/*
 * Copyright 2013-2023, Seqera Labs
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
import groovy.util.logging.Slf4j

/**
 * Models container engine configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
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

    ContainerConfig setEnabled(boolean value) {
        put('enabled', value)
        return this
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

    String fusionOptions() {
        final result = get('fusionOptions')
        return result!=null ? result : defaultFusionOptions()
    }

    protected String defaultFusionOptions() {
        final eng = getEngine()
        if( !eng )
            return null
        if( eng=='docker' || eng=='podman' )
            return '--rm --privileged'
        if( eng=='singularity' || eng=='apptainer' )
            return null
        log.warn "Fusion file system is not supported by '$eng' container engine"
        return null
    }
}
