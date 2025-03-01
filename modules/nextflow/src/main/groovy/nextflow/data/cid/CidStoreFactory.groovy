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
 *
 */

package nextflow.data.cid

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.data.config.DataConfig
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Factory for CidStore
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
abstract class CidStoreFactory implements ExtensionPoint {

    protected abstract CidStore newInstance(DataConfig config)

    static CidStore create(DataConfig config){
        final all = Plugins.getPriorityExtensions(CidStoreFactory)
        if( !all )
            throw new IllegalStateException("Unable to find Nextflow CID store factory")
        final factory = all.first()
        log.debug "Using Nextflow CID store factory: ${factory.getClass().getName()}"
        return factory.newInstance(config)


    }


}
