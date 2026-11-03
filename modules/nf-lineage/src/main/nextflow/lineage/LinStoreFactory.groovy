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

package nextflow.lineage

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.lineage.config.LineageConfig
import nextflow.plugin.Plugins
import nextflow.util.TestOnly
import org.pf4j.ExtensionPoint

/**
 * Factory for {@link LinStore} objects
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
abstract class LinStoreFactory implements ExtensionPoint {

    private static LinStore instance

    private static boolean initialized

    protected abstract boolean canOpen(LineageConfig config)

    protected abstract LinStore newInstance(LineageConfig config)

     static LinStore create(LineageConfig config){
         final factory = Plugins
             .getPriorityExtensions(LinStoreFactory)
             .find( f-> f.canOpen(config))
         if( !factory )
            throw new IllegalStateException("Unable to find Nextflow Lineage store factory")
        log.debug "Using Nextflow Lineage store factory: ${factory.getClass().getName()}"
        return factory.newInstance(config)
    }

    static LinStore getOrCreate(Session session) {
        if( instance || initialized )
            return instance
        synchronized (LinStoreFactory.class) {
            if( instance || initialized )
                return instance
            initialized = true
            final config = LineageConfig.create(session)
            if( !config.enabled )
                return null
            return instance = create(config)
        }
    }

    @TestOnly
    static void reset(){
        synchronized (LinStoreFactory.class) {
            instance = null
            initialized = false
        }
    }

}
