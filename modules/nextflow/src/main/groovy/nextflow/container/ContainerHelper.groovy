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
package nextflow.container

import groovy.transform.CompileStatic
import nextflow.SysEnv

@CompileStatic
class ContainerHelper {

    static boolean entrypointOverride() {
        return SysEnv.getBool('NXF_CONTAINER_ENTRYPOINT_OVERRIDE', false)
    }

    static boolean fixOwnership(ContainerConfig config) {
        return config instanceof DockerConfig && config.fixOwnership
    }

    static List<String> parseEnvWhitelist(Object value) {
        if( !value )
            return []
        if( value instanceof CharSequence )
            return value.tokenize(',').collect { it.trim() }
        if( value instanceof List )
            return value as List<String>
        throw new IllegalArgumentException("Not a valid `envWhitelist` argument: $value")
    }

}
