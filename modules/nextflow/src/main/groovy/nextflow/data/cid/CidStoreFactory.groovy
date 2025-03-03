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
import nextflow.Session
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

    private static CidStore instance

    private static boolean initialized

    protected abstract CidStore newInstance(DataConfig config)

    private static CidStore create(DataConfig config){
        final all = Plugins.getPriorityExtensions(CidStoreFactory)
        if( !all )
            throw new IllegalStateException("Unable to find Nextflow CID store factory")
        final factory = all.first()
        log.debug "Using Nextflow CID store factory: ${factory.getClass().getName()}"
        return factory.newInstance(config)
    }

    static CidStore getOrCreate(Session session) {
        if( instance || initialized )
            return instance
        synchronized (session) {
            if( instance || initialized )
                return instance
            initialized = true
            final config = DataConfig.create(session)
            if( !config.enabled )
                return null
            return instance = create(config)
        }
    }

}
