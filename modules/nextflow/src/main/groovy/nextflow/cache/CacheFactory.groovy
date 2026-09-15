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

package nextflow.cache

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Factory class that create an instance of the {@link CacheDB}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class CacheFactory implements ExtensionPoint {

    protected abstract CacheDB newInstance(Session session, UUID uniqueId, String runName, Path home=null)

    /**
     * Whether this factory applies to the given session.
     *
     * @param session The current {@link Session}, or null when running outside of a workflow
     *        execution e.g. the `log` and `clean` commands
     */
    protected boolean enabled(Session session) { true }

    static CacheDB create(Session session, UUID uniqueId, String runName, Path home=null) {
        final all = Plugins.getPriorityExtensions(CacheFactory)
        final factory = all.find(it -> it.enabled(session))
        if( !factory )
            throw new IllegalStateException("Unable to find Nextflow cache factory")
        log.debug "Using Nextflow cache factory: ${factory.getClass().getName()}"
        return factory.newInstance(session, uniqueId, runName, home)
    }

}
