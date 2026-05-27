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

package nextflow.module

import groovy.transform.CompileStatic
import io.seqera.http.HxClient
import io.seqera.npr.client.RegistryClient
import nextflow.config.RegistryConfig
import nextflow.util.RetryConfig
import nextflow.util.TestOnly

import java.net.http.HttpClient

/**
 * Factory for creating a {@link RegistryClient} configured from Nextflow's {@link RegistryConfig}.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@CompileStatic
class RegistryClientFactory {

    private static volatile RegistryClient instance

    /** Default registry always appended as a fallback for module resolution */
    private static String defaultRegistry = RegistryConfig.DEFAULT_REGISTRY_URL

    @TestOnly
    static void setDefaultRegistry(String url) {
        defaultRegistry = url
    }

    static RegistryClient forConfig(RegistryConfig config) {
        final cfg = config ?: new RegistryConfig()
        // the default registry is always queried first as the baseline,
        // followed by any registries configured in the `registry` scope
        final urls = new ArrayList<String>()
        if( defaultRegistry )
            urls.add(defaultRegistry)
        for( String it : cfg.allUrls )
            if( !urls.contains(it) )
                urls.add(it)
        return new RegistryClient(
            urls,
            cfg.apiKey,
            HxClient.newBuilder()
                .retryConfig(RetryConfig.config())
                .followRedirects(HttpClient.Redirect.NORMAL)
                .build()
        )
    }

    static RegistryClient instance(RegistryConfig config) {
        if( instance == null ) {
            synchronized (RegistryClientFactory) {
                if( instance == null )
                    instance = forConfig(config)
            }
        }
        return instance
    }

    static void reset() {
        instance = null
    }
}
