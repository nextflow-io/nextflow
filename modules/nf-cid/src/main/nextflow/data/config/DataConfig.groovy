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

package nextflow.data.config

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session

/**
 * Model workflow data config
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataConfig {

    final DataStoreOpts store

    final boolean enabled

    DataConfig(Map opts) {
        this.store = new DataStoreOpts(opts.store as Map ?: Map.of())
        this.enabled = opts.enabled as boolean ?: false
    }

    static Map<String,Object> asMap() {
        session?.config?.navigate('workflow.data') as Map ?: new HashMap<String,Object>()
    }

    static DataConfig create(Session session) {
        if( session ) {
            return new DataConfig( session.config.navigate('workflow.data') as Map ?: Map.of())
        }
        else
            throw new IllegalStateException("Missing Nextflow session")
    }

    static DataConfig create() {
        create(getSession())
    }

    static private Session getSession() {
        return Global.session as Session
    }
}
