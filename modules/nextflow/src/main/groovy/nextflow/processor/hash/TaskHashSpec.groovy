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
import nextflow.util.CacheHelper
import nextflow.util.EncodingRules

/**
 * A task hash version expressed as data: which keys are hashed, in what order, how
 * each is extracted, and how values are encoded.
 *
 * A spec is immutable and, once shipped, frozen: editing one changes the hashes of
 * every task already recorded under it. A new behaviour is a new spec.
 */
@CompileStatic
class TaskHashSpec {

    final String id

    final Map<HashKey,Contributor> bindings

    final EncodingRules encoding

    private volatile String fingerprintValue

    TaskHashSpec(String id, Map<HashKey,Contributor> bindings, EncodingRules encoding) {
        this.id = id
        this.bindings = Collections.unmodifiableMap(new LinkedHashMap<HashKey,Contributor>(bindings))
        this.encoding = encoding
    }

    List<HashKey> keys() {
        return new ArrayList<HashKey>(bindings.keySet())
    }

    /**
     * Stable text form from which the fingerprint is derived. It names only semantics
     * — id, ordered keys, extraction identity, encoding, hash function — so that
     * refactoring implementation detail cannot move a historical spec's identity.
     * Never reformat this method.
     */
    String canonicalForm() {
        final sb = new StringBuilder()
        sb.append('id=').append(id).append('\n')
        for( Map.Entry<HashKey,Contributor> b : bindings.entrySet() ) {
            sb.append('key=').append(b.key.name()).append(':').append(b.value.canonicalName()).append('\n')
        }
        sb.append('encoding=').append(encoding.canonicalForm()).append('\n')
        sb.append('function=murmur3_128').append('\n')
        return sb.toString()
    }

    String fingerprint() {
        if( fingerprintValue == null ) {
            fingerprintValue = CacheHelper.hasher(canonicalForm()).hash().toString()
        }
        return fingerprintValue
    }

    @Override
    String toString() {
        return "TaskHashSpec[${id}@${fingerprint()}]"
    }

    static final String SUPPORTED_FUNCTION = 'murmur3_128'

    static TaskHashSpec fromJson(String json, String origin) {
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
                boolAt(encodingNode as Map, 'cacheFunnelFirst', id),
                boolAt(encodingNode as Map, 'assetRootDetection', id))

        final keysNode = root.get('keys')
        if( !(keysNode instanceof List) || !keysNode ) {
            throw new IllegalArgumentException("Task hash spec ${id} is missing 'keys'")
        }

        final bindings = new LinkedHashMap<HashKey,Contributor>()
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
            if( bindings.containsKey(key) ) {
                throw new IllegalArgumentException("Task hash spec ${id} binds ${keyName} more than once")
            }
            bindings.put(key, Contributors.get(contributor))
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
