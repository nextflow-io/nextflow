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

package nextflow.trace

import spock.lang.Specification

class RuntimeMetadataStoreTest extends Specification {

    def 'should put and get values' () {
        given:
        def store = new RuntimeMetadataStore()

        when:
        store.put('foo', 'bar')
        store.put('count', 42)

        then:
        store.get('foo') == 'bar'
        store.get('count') == 42
        store.get('missing') == null
    }

    def 'should remove key when value is null' () {
        given:
        def store = new RuntimeMetadataStore()
        store.put('foo', 'bar')

        when:
        store.put('foo', null)

        then:
        store.get('foo') == null
        store.isEmpty()
    }

    def 'should overwrite existing values' () {
        given:
        def store = new RuntimeMetadataStore()

        when:
        store.put('foo', 'one')
        store.put('foo', 'two')

        then:
        store.get('foo') == 'two'
    }

    def 'should reject null or empty key' () {
        given:
        def store = new RuntimeMetadataStore()

        when:
        store.put(null, 'x')
        then:
        thrown(IllegalArgumentException)

        when:
        store.put('', 'x')
        then:
        thrown(IllegalArgumentException)
    }

    def 'snapshot should be immutable' () {
        given:
        def store = new RuntimeMetadataStore()
        store.put('foo', 'bar')

        when:
        def snap = store.snapshot()
        snap.put('x', 'y')

        then:
        thrown(UnsupportedOperationException)
    }

    def 'snapshot should not reflect later mutations' () {
        given:
        def store = new RuntimeMetadataStore()
        store.put('foo', 'bar')

        when:
        def snap = store.snapshot()
        store.put('foo', 'changed')
        store.put('new', 'value')

        then:
        snap == [foo: 'bar']
    }

    def 'isEmpty should reflect store state' () {
        given:
        def store = new RuntimeMetadataStore()

        expect:
        store.isEmpty()

        when:
        store.put('foo', 'bar')
        then:
        !store.isEmpty()

        when:
        store.put('foo', null)
        then:
        store.isEmpty()
    }
}
