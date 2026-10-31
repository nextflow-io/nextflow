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

package nextflow.config;

import java.util.ServiceLoader;

/**
 * Creates the Nextflow config parser, which lives in the runtime and is registered as a service
 */
public interface ConfigParserProvider {

    ConfigParser create();

    static ConfigParser createParser() {
        // this module's class loader, not the thread context one, so plugin threads find the service too
        final ConfigParserProvider provider = ServiceLoader.load(ConfigParserProvider.class, ConfigParserProvider.class.getClassLoader()).findFirst().orElse(null);
        if( provider == null )
            throw new IllegalStateException("No Nextflow config parser available -- the nextflow runtime must be on the classpath");
        return provider.create();
    }
}
