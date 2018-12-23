/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.processor

import groovy.transform.CompileStatic

/**
 * Context shared across multiple {@link TaskHandler} instances of the same type.
 *
 * This is required to aggregate multiple remote API invocations and execute them
 * as a sole remote request.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class BatchContext<K,V> {

    private List<K> entries = new ArrayList()

    private Map<K, V> cache = new HashMap<>()

    int size() { entries.size() }

    void collect( K key ) {
        assert key != null
        if( !entries.contains(key) )
            entries.add(key)
    }

    Collection<K> getBatchFor( K key, int n ) {
        def p = entries.indexOf(key)
        if( p == -1 )
            throw new IllegalStateException("Cannot find key `$key` in batch collector")
        entries.subList(p, Math.min(p+n, entries.size()))
    }

    def V get( K key ) {
        (V)cache.get(key)
    }

    boolean contains( K key ) {
        cache.get(key) != null
    }

    void put( K key, V payload ) {
        cache.put(key, payload)
    }

}
