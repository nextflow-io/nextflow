/*
 * Copyright 2021, Sage-Bionetworks
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
 *
 */

package nextflow.secret

import groovy.util.logging.Slf4j
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins

/**
 * Implements dynamic secret providing loading strategy
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@Singleton
class SecretsLoader {

    private SecretsProvider provider

    static boolean isEnabled() {
        SysEnv.get('NXF_ENABLE_SECRETS', 'true') == 'true'
    }

    SecretsProvider load() {
        if( provider )
            return provider
        synchronized (this) {
            if( provider )
                return provider
            return provider = load0()
        }
    }

    private SecretsProvider load0() {
        // discover all available secrets provider
        final all = Plugins.getPriorityExtensions(SecretsProvider)
        // find first activable in the current environment
        final provider = all.find { it.activable() }
        log.debug "Discovered secrets providers: $all - activable => $provider"
        if( provider )
            return provider.load()
        else
            throw new AbortOperationException("Unable to load secrets provider")
    }

    void reset() {
        provider=null
    }

    static protected makeSecretsContext(SecretsProvider provider) {

        return new Object() {
            def getProperty(String name) {
                if( !provider )
                    throw new AbortOperationException("Unable to resolve secrets.$name - no secret provider is available")
                provider.getSecret(name)?.value
            }
        }
    }

    static Object secretContext() {
        final provider = isEnabled() ? getInstance().load() : new NullProvider()
        return makeSecretsContext(provider)
    }
    
}
