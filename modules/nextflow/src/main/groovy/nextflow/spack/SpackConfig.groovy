/*
 * Copyright 2022-2023, Pawsey Supercomputing Research Centre
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

package nextflow.spack

import groovy.transform.CompileStatic

/**
 * Model Spack configuration
 * 
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@CompileStatic
class SpackConfig extends LinkedHashMap {

    private Map<String,String> env

    /* required by Kryo deserialization -- do not remove */
    private SpackConfig() { }

    SpackConfig(Map config, Map<String, String> env) {
        super(config)
        this.env = env
    }

    boolean isEnabled() {
        def enabled = get('enabled')
        if( enabled == null )
            enabled = env.get('NXF_SPACK_ENABLED')
        return enabled?.toString() == 'true'
    }

    List<String> getChannels() {
        final value = get('channels')
        if( !value ) {
            return Collections.<String>emptyList()
        }
        if( value instanceof List ) {
            return value
        }
        if( value instanceof CharSequence ) {
            return value.tokenize(',').collect(it -> it.trim())
        }

        throw new IllegalArgumentException("Unexpected spack.channels value: $value")
    }
}
