/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
