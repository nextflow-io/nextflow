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

package nextflow.config

import groovy.transform.CompileStatic
import nextflow.secret.SecretHolder
import nextflow.secret.SecretsProvider

/**
 * Represent Nextflow config as Map
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConfigMap extends LinkedHashMap {

    ConfigMap() {
    }

    ConfigMap(int initialCapacity) {
        super(initialCapacity)
    }

    ConfigMap(Map content) {
        super(content)
    }

    @Override
    Object get(Object key) {
        final result = super.get(key)
        // check if it's a secret value
        if( result instanceof SecretHolder && result.isBound() ) {
            return result.call()
        }
        return result
    }

    ConfigMap withSecretProvider(SecretsProvider provider) {
        withSecretProvider0(provider,this)
        return this
    }

    private withSecretProvider0(SecretsProvider provider, Map map) {
        for( Object key : map.keySet() ) {
            def entry = map.get(key)
            // traverse nested config map objects
            if( entry instanceof Map ) {
                withSecretProvider0(provider, entry)
            }
            // look for all secret holders in the config map
            // and bind the secrets provider
            if( entry instanceof SecretHolder ) {
                entry.bind(provider)
            }
            // same bind secret holders in Gstring objects
            else if( entry instanceof GString ) {
                final str = (GString)entry
                for( Object value : str.getValues() ) {
                    if( value instanceof SecretHolder ) {
                        value.bind(provider)
                    }
                }
            }
        }
    }
}
