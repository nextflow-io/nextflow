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

package nextflow.plugin

import groovy.transform.CompileStatic
import org.pf4j.update.UpdateRepository

/**
 * Extension to pf4j's UpdateRepository which supports pre-fetching
 * metadata for a specified set of plugins.
 *
 * This gives the ability to avoid downloading metadata for unused
 * plugins.
 */
@CompileStatic
interface PrefetchUpdateRepository extends UpdateRepository {
    /**
     * This will be called when Nextflow starts, before
     * initialising the plugins.
     */
    void prefetch(List<PluginRef> plugins)
}
