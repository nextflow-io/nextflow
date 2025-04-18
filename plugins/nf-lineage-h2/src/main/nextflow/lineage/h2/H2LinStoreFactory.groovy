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
 *
 */

package nextflow.lineage.h2

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.config.LineageConfig
import nextflow.plugin.Priority

@Slf4j
@CompileStatic
@Priority(-10)  // <-- lower is higher, this is needed to override default provider behavior
class H2LinStoreFactory extends LinStoreFactory {

    @Override
    boolean canOpen(LineageConfig config) {
        return config.store.location.startsWith('jdbc:h2:')
    }

    @Override
    protected LinStore newInstance(LineageConfig config) {
        return new H2LinStore().open(config)
    }
}
