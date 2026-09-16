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

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
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

    final List<KeyBinding> bindings

    final EncodingRules encoding

    private volatile String fingerprintValue

    TaskHashSpec(String id, List<KeyBinding> bindings, EncodingRules encoding) {
        this.id = id
        this.bindings = Collections.unmodifiableList(new ArrayList<KeyBinding>(bindings))
        this.encoding = encoding
    }

    List<HashKey> keys() {
        return bindings.collect { KeyBinding it -> it.key }
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
        for( KeyBinding b : bindings ) {
            sb.append('key=').append(b.key.name()).append(':').append(b.contributor.canonicalName()).append('\n')
        }
        sb.append('encoding=').append(encoding.canonicalForm()).append('\n')
        sb.append('function=murmur3_128').append('\n')
        return sb.toString()
    }

    String fingerprint() {
        if( fingerprintValue == null ) {
            fingerprintValue = Hashing.murmur3_128().newHasher()
                .putUnencodedChars(canonicalForm())
                .hash()
                .toString()
        }
        return fingerprintValue
    }

    @Override
    String toString() {
        return "TaskHashSpec[${id}@${fingerprint()}]"
    }
}
