/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

    /* required by Kryo deserialization -- do not remove */
    private ContainerConfig() { }

    ContainerConfig(Map config) {
        super(config)
    }

    boolean isEnabled() {
        get('enabled')?.toString() == 'true'
    }

    boolean isLegacy() {
        get(legacy)?.toString() == 'true'
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
}
