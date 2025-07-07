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

package nextflow.plugin

import com.google.common.hash.Hasher
import groovy.transform.Canonical
import nextflow.util.CacheFunnel
import nextflow.util.CacheHelper
/**
 * Model a plugin Id and version
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
class PluginSpec implements CacheFunnel, Comparable<PluginSpec> {

    /**
     * Plugin unique ID
     */
    String id

    /**
     * The plugin version
     */
    String version

    /**
     * Parse a plugin fully-qualified ID eg. nf-amazon@1.2.0
     *
     * @param fqid The fully qualified plugin id
     * @return A {@link PluginSpec} representing the plugin
     */
    static PluginSpec parse(String fqid, DefaultPlugins defaultPlugins=null) {
        final tokens = fqid.tokenize('@') as List<String>
        final id = tokens[0]
        final ver = tokens[1]
        if( ver || defaultPlugins==null )
            return new PluginSpec(id, ver)
        if( defaultPlugins.hasPlugin(id) )
            return defaultPlugins.getPlugin(id)
        return new PluginSpec(id)
    }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(id)
        if( version )
            hasher.putUnencodedChars(version)
        return hasher
    }

    @Override
    String toString() {
        version ? "${id}@${version}" : id
    }

    @Override
    int compareTo(PluginSpec that) {
        return this.toString() <=> that.toString()
    }
}
