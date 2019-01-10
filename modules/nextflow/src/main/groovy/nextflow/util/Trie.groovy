/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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



