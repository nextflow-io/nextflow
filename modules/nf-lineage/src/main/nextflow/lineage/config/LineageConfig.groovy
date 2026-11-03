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

package nextflow.lineage.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.Global
import nextflow.Session
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Model workflow data lineage config
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("lineage")
@Description("""
    The `lineage` scope controls the generation of lineage metadata.
""")
@ToString
@CompileStatic
class LineageConfig implements ConfigScope {

    final LineageStoreOpts store

    @ConfigOption
    @Description("""
        Enable generation of lineage metadata (default: `false`).
    """)
    final boolean enabled

    /* required by extension point -- do not remove */
    LineageConfig() {}

    LineageConfig(Map opts) {
        this.store = new LineageStoreOpts(opts.store as Map ?: Collections.emptyMap())
        this.enabled = opts.enabled as boolean
    }

    static Map<String,Object> asMap() {
        final result = session?.config?.navigate('lineage') as Map
        return result != null ? result : new HashMap<String,Object>()
    }

    static LineageConfig create(Map config) {
        new LineageConfig(config.lineage as Map ?: Collections.emptyMap())
    }

    static LineageConfig create(Session session) {
        if( session ) {
            return LineageConfig.create(session.config)
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
