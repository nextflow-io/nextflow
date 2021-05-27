/*
 * Copyright 2020-2021, Seqera Labs
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
class PluginSpec implements CacheFunnel {

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
    static PluginSpec parse(String fqid) {
        final tokens = fqid.tokenize('@') as List<String>
        final id = tokens[0]
        return new PluginSpec(id, tokens[1])
    }

    @Override
    Hasher funnel(Hasher hasher, CacheHelper.HashMode mode) {
        hasher.putUnencodedChars(id)
        if( version )
            hasher.putUnencodedChars(version)
        return hasher
    }

    String toString() {
        version ? "${id}@${version}" : id
    }
}
