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

package nextflow.cache

import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
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

    protected abstract CacheDB newInstance(UUID uniqueId, String runName, Path home=null)

    static CacheDB create(UUID uniqueId, String runName, Path home=null) {
        final all = Plugins.getPriorityExtensions(CacheFactory)
        if( !all )
            throw new IllegalStateException("Unable to find Nextflow cache factory")
        final factory = all.first()
        log.debug "Using Nextflow cache factory: ${factory.getClass().getName()}"
        return factory.newInstance(uniqueId, runName, home)
    }

}
