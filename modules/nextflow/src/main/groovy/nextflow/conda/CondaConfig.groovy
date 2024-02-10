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

package nextflow.conda

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 * Model Conda configuration
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CondaConfig extends LinkedHashMap {

    private Map<String,String> env

    /* required by Kryo deserialization -- do not remove */
    private CondaConfig() { }

    CondaConfig(Map config, Map<String, String> env) {
        super(config)
        this.env = env
    }

    boolean isEnabled() {
        def enabled = get('enabled')
        if( enabled == null )
            enabled = env.get('NXF_CONDA_ENABLED')
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

        throw new IllegalArgumentException("Unexpected conda.channels value: $value")
    }

    Duration createTimeout() {
        get('createTimeout') as Duration
    }

    String createOptions() {
        get('createOptions') as String
    }

    Path cacheDir() {
        get('cacheDir') as Path
    }

    boolean useMamba() {
        get('useMamba') as boolean
    }

    boolean useMicromamba() {
        get('useMicromamba') as boolean
    }
}
