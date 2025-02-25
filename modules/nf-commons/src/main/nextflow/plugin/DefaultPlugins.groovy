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

import groovy.transform.CompileStatic

/**
 * Model the collection of default plugins used if no plugins are provided by the user
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DefaultPlugins {

    public static final DefaultPlugins INSTANCE = new DefaultPlugins()

    private Map<String,PluginSpec> plugins = new HashMap<>(20)

    protected DefaultPlugins() {
        final meta = this.class.getResourceAsStream('/META-INF/plugins-info.txt')?.text
        plugins = parseMeta(meta)
    }

    protected Map<String,PluginSpec> parseMeta(String meta) {
        if( !meta )
            return Collections.emptyMap()

        final result = new HashMap(20)
        for( String line : meta.readLines() ) {
            final spec = PluginSpec.parse(line)
            result[spec.id] = spec
        }
        return result
    }

    PluginSpec getPlugin(String pluginId) throws IllegalArgumentException {
        if( !pluginId )
            throw new IllegalArgumentException("Missing pluginId argument")
        final result = plugins.get(pluginId)
        if( !result )
            throw new IllegalArgumentException("Unknown Nextflow plugin '$pluginId'")
        return result
    }

    boolean hasPlugin(String pluginId) {
        return plugins.containsKey(pluginId)
    }

    List<PluginSpec> getPlugins() {
        return new ArrayList(plugins.values())
    }

    String toSortedString(String divisor=',') {
        getPlugins().sort().join(divisor)
    }

    @Override
    String toString() {
        return "DefaultPlugins${getPlugins()}"
    }
}
