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

import nextflow.util.EncodingRules
import spock.lang.Specification

class TaskHashSpecTest extends Specification {

    static Contributor contrib(String name) {
        return new Contributor() {
            @Override
            String canonicalName() { return name }

            @Override
            List<Object> emit(HashContext ctx) { return [] }
        }
    }

    def 'canonical form names the id, the ordered keys with their extraction, and the encoding'() {
        given:
        def spec = new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.SESSION_ID, contrib('sessionId')),
            new KeyBinding(HashKey.TASK_SOURCE, contrib('taskSource'))
        ], EncodingRules.RECORD_TYPES)

        expect:
        spec.canonicalForm() == '''\
id=std/test
key=SESSION_ID:sessionId
key=TASK_SOURCE:taskSource
encoding=orderIndependentMaps=true;cacheFunnelFirst=true
function=murmur3_128
'''
    }

    def 'fingerprint is stable, and moves when any part of the canonical form moves'() {
        given:
        def base = new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.RECORD_TYPES)

        expect: 'stable across calls and across equal specs'
        base.fingerprint() == base.fingerprint()
        and:
        base.fingerprint() == new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.RECORD_TYPES).fingerprint()

        and: 'a different extraction for the same key moves it'
        base.fingerprint() != new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId.other'))], EncodingRules.RECORD_TYPES).fingerprint()

        and: 'different encoding moves it'
        base.fingerprint() != new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.LEGACY).fingerprint()

        and: 'different key order moves it'
        new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.SESSION_ID, contrib('a')),
            new KeyBinding(HashKey.TASK_SOURCE, contrib('b'))
        ], EncodingRules.RECORD_TYPES).fingerprint() != new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.TASK_SOURCE, contrib('b')),
            new KeyBinding(HashKey.SESSION_ID, contrib('a'))
        ], EncodingRules.RECORD_TYPES).fingerprint()
    }

    def 'keys() reports the declared keys in order'() {
        given:
        def spec = new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.TASK_SOURCE, contrib('a')),
            new KeyBinding(HashKey.CONDA, contrib('b'))
        ], EncodingRules.RECORD_TYPES)

        expect:
        spec.keys() == [HashKey.TASK_SOURCE, HashKey.CONDA]
    }
}
