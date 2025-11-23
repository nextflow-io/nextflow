/*
 * Copyright 2013-2025, Seqera Labs
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

package io.seqera.tower.plugin.scm

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.scm.jgit.TransportSeqera
import nextflow.plugin.Priority
import nextflow.scm.GitUrl
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryFactory
import nextflow.scm.RepositoryProvider

import java.util.concurrent.atomic.AtomicBoolean

/**
 * Implements a factory to create an instance of {@link SeqeraRepositoryProvider}
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@Priority(-10)
@CompileStatic
class SeqeraRepositoryFactory extends RepositoryFactory {

    private static AtomicBoolean registered = new AtomicBoolean(false)

    @Override
    protected RepositoryProvider createProviderInstance(ProviderConfig config, String project) {
        if (!registered.get()) {
            registered.set(true)
            TransportSeqera.register()
        }

        return config.platform == 'seqera'
            ? new SeqeraRepositoryProvider(project, config)
            : null
    }

    @Override
    protected ProviderConfig getConfig(List<ProviderConfig> providers, GitUrl url) {
        // Only handle seqera:// URLs
        if (url.protocol != 'seqera') {
            return null
        }

        // Seqera repository config depends on the data-link ID stored as domain
        def config = providers.find { it -> it.domain == url.domain }
        if (config) {
            log.debug "Git url=$url (1) -> config=$config"
            return config
        }

        // Create a new instance if not found
        config = new SeqeraProviderConfig(url.domain)

        return config
    }

    @Override
    protected ProviderConfig createConfigInstance(String name, Map attrs) {
        final copy = new HashMap(attrs)
        return copy.platform == 'seqera'
            ? new SeqeraProviderConfig(name, copy)
            : null
    }
}
