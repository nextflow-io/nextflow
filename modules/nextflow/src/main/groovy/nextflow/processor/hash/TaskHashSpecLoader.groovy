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
package nextflow.processor.hash

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import nextflow.util.EncodingRules

/**
 * Builds a {@link TaskHashSpec} from its JSON description.
 *
 * Every failure here is fatal by design. A spec file names the keys whose values become
 * a cache key, so a typo that silently resolved to something else — or to nothing —
 * would produce hashes that cannot be accounted for afterwards. There is no lenient mode.
 */
@CompileStatic
class TaskHashSpecLoader {

    static final String SUPPORTED_FUNCTION = 'murmur3_128'

    static TaskHashSpec load(String json, String origin) {
        final parsed = new JsonSlurper().parseText(json)
        if( !(parsed instanceof Map) ) {
            throw new IllegalArgumentException("Task hash spec ${origin} must be a JSON object")
        }
        final root = parsed as Map

        final id = root.get('id') as String
        if( !id ) {
            throw new IllegalArgumentException("Task hash spec ${origin} is missing 'id'")
        }

        final function = root.get('function') as String
        if( function != SUPPORTED_FUNCTION ) {
            throw new IllegalArgumentException("Task hash spec ${id} declares unsupported hash function '${function}' -- only '${SUPPORTED_FUNCTION}' is implemented")
        }

        final encodingNode = root.get('encoding')
        if( !(encodingNode instanceof Map) ) {
            throw new IllegalArgumentException("Task hash spec ${id} is missing 'encoding'")
        }
        final encoding = new EncodingRules(
                boolAt(encodingNode as Map, 'orderIndependentMaps', id),
                boolAt(encodingNode as Map, 'cacheFunnelFirst', id))

        final keysNode = root.get('keys')
        if( !(keysNode instanceof List) || !keysNode ) {
            throw new IllegalArgumentException("Task hash spec ${id} is missing 'keys'")
        }

        final bindings = new ArrayList<KeyBinding>()
        final seen = new HashSet<HashKey>()
        for( Object entry : (keysNode as List) ) {
            if( !(entry instanceof Map) ) {
                throw new IllegalArgumentException("Task hash spec ${id} has a malformed entry in 'keys'")
            }
            final e = entry as Map
            final keyName = e.get('key') as String
            final contributor = e.get('contributor') as String
            if( !keyName || !contributor ) {
                throw new IllegalArgumentException("Task hash spec ${id} has an entry missing 'key' or 'contributor'")
            }
            final key = parseKey(keyName, id)
            if( !seen.add(key) ) {
                throw new IllegalArgumentException("Task hash spec ${id} binds ${keyName} more than once")
            }
            bindings.add(new KeyBinding(key, ContributorRegistry.get(contributor)))
        }

        return new TaskHashSpec(id, bindings, encoding)
    }

    private static HashKey parseKey(String name, String specId) {
        try {
            return HashKey.valueOf(name)
        }
        catch( IllegalArgumentException e ) {
            throw new IllegalArgumentException("Task hash spec ${specId} names unknown hash key '${name}' -- available: ${HashKey.values()*.name().join(', ')}", e)
        }
    }

    private static boolean boolAt(Map node, String field, String specId) {
        final value = node.get(field)
        if( !(value instanceof Boolean) ) {
            throw new IllegalArgumentException("Task hash spec ${specId} encoding field '${field}' must be true or false")
        }
        return (Boolean) value
    }
}
