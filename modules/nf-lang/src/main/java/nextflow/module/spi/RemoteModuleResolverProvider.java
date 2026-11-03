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

package nextflow.module.spi;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.ServiceLoader;

/**
 * Provider for accessing RemoteModuleResolver implementations via SPI.
 *
 * <p>This class uses the Java ServiceLoader mechanism to discover and load
 * implementations of RemoteModuleResolver. It selects the implementation
 * with the highest priority.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class RemoteModuleResolverProvider {

    private static final Logger log = LoggerFactory.getLogger(RemoteModuleResolverProvider.class);
    private static RemoteModuleResolver instance;

    /**
     * Get the RemoteModuleResolver instance with the highest priority.
     *
     * <p>This method lazily loads and caches the resolver. It discovers all
     * implementations via ServiceLoader and selects the one with the highest
     * priority value.
     *
     * <p>If no implementations are found, returns the FallbackRemoteModuleResolver
     * which throws an informative exception.
     *
     * @return The RemoteModuleResolver instance with highest priority
     */
    public static synchronized RemoteModuleResolver getInstance() {
        if (instance == null) {
            instance = loadResolver();
        }
        return instance;
    }

    private static RemoteModuleResolver loadResolver() {
        List<RemoteModuleResolver> resolvers = new ArrayList<>();
        ServiceLoader<RemoteModuleResolver> loader = ServiceLoader.load(RemoteModuleResolver.class);

        // Collect all available resolvers
        for (RemoteModuleResolver resolver : loader) {
            resolvers.add(resolver);
            log.debug("Discovered RemoteModuleResolver: {} with priority {}",
                resolver.getClass().getName(), resolver.getPriority());
        }

        // Sort by priority (highest first)
        resolvers.sort(Comparator.comparingInt(RemoteModuleResolver::getPriority).reversed());

        if (resolvers.isEmpty()) {
            log.warn("No RemoteModuleResolver implementations found via SPI, using fallback");
            return new FallbackRemoteModuleResolver();
        }

        RemoteModuleResolver selected = resolvers.get(0);
        log.debug("Selected RemoteModuleResolver: {} with priority {}",
            selected.getClass().getName(), selected.getPriority());

        return selected;
    }

    /**
     * Reset the cached instance. Used primarily for testing.
     */
    public static synchronized void reset() {
        instance = null;
    }
}