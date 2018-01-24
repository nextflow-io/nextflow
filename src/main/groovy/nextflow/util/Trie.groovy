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

package nextflow.util

import groovy.transform.CompileStatic


/**
 * Implements a basic Trie data structure
 *
 * @link http://en.wikipedia.org/wiki/Trie
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Trie<T> {

    private T vertex

    private List<Trie<T>> children

    T getVertex() { vertex }

    List<Trie<T>> getChildren() { children }

    Trie( T obj ) {
        this.vertex = obj
    }

    Trie append( T value ) {

        if( children == null ) {
            def result = new Trie<T>(value)
            children = new LinkedList<Trie<T>>()
            children << result
            return result
        }

        def result = children?.find { node -> node.vertex == value }
        if( !result ) {
            result = new Trie<T>(value)
            children << result
            return result
        }

        return result
    }

    Trie append( List<T> values ) {
        if(!values)
            return null

        def v = values.head()
        def node = append(v)
        node.append( values.tail() )
        return node
    }

    Trie append( T... values ) {
        append( values as List<T> )
    }

    Trie getChild( T value ) {
        children?.find { node -> node.vertex == value }
    }

    List<T> longest() {
        longestImpl new LinkedList<T>()
    }

    private List<T> longestImpl(List<T> result) {
        result << vertex
        if( children?.size() == 1 ) {
            children.get(0).longestImpl(result)
        }
        return result
    }



}



