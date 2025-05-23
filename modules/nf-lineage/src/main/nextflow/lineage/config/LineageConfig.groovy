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

package nextflow.lineage.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session

/**
 * Model workflow data lineage config
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString
@CompileStatic
class LineageConfig {

    final LineageStoreOpts store

    final boolean enabled

    LineageConfig(Map opts) {
        this.store = new LineageStoreOpts(opts.store as Map ?: Map.of())
        this.enabled = opts.enabled as boolean ?: false
    }

    static Map<String,Object> asMap() {
        final result = session?.config?.navigate('lineage') as Map
        return result != null ? result : new HashMap<String,Object>()
    }

    static LineageConfig create(Session session) {
        if( session ) {
            return new LineageConfig( session.config.navigate('lineage') as Map ?: Map.of())
        }
        else
            throw new IllegalStateException("Missing Nextflow session")
    }

    static LineageConfig create() {
        create(getSession())
    }

    static private Session getSession() {
        return Global.session as Session
    }
}
