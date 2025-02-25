/*
 * Copyright 2013-2024, Seqera Labs
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

    private T node

    private List<Trie<T>> children

    T getNode() { node }

    List<Trie<T>> getChildren() { children }

    Trie( T obj ) {
        this.node = obj
    }

    protected Trie addNode( T node ) {

        if( children == null ) {
            final result = new Trie<T>(node)
            children = new LinkedList<Trie<T>>()
            children.add(result)
            return result
        }

        def result = children?.find { trie -> trie.node == node }
        if( !result ) {
            result = new Trie<T>(node)
            children.add(result)
            return result
        }

        return result
    }

    Trie addPath( List<T> nodes ) {
        if(!nodes)
            return null

        def node = addNode(nodes.head())
        node.addPath(nodes.tail())
        return node
    }

    Trie addPath( T... nodes ) {
        addPath( nodes as List<T> )
    }

    Trie getChild( T name ) {
        children?.find { trie -> trie.node == name }
    }

    List<T> longest() {
        longestImpl new LinkedList<T>()
    }

    private List<T> longestImpl(List<T> result) {
        result << node
        if( children?.size() == 1 ) {
            children.get(0).longestImpl(result)
        }
        return result
    }

    List<List<T>> traverse(T stop=null) {
        def result = new ArrayList<List<T>>()
        traverse0(new LinkedList<T>(), result, stop)
        return result
    }

    private void traverse0(List<T> current, List<List<T>> result, T stop) {
        current.add(node)
        if( !children || children.any { it.node==stop } ) {
            result.add(current)
            return
        }
        for( Trie t : children ) {
            t.traverse0(new ArrayList<T>(current), result, stop)
        }
    }

}



