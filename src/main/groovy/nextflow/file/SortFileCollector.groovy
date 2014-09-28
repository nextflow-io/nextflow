/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.transform.ToString
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class SortFileCollector<K, V> implements Map<K,V>, Iterable<K> {

    private int count

    private Closure<Comparable> comparator

    private Map<Integer,V> store = new HashMap<>()

    private Map<Integer,BTree<K>> index = new HashMap<>()

    @PackageScope
    Map<Integer,V> getStore() { store }

    @PackageScope
    Map<Integer,BTree<K>> getIndex() { index }

    @PackageScope
    int getCount() { count }

    protected int compare( K k1, K k2 ) {
        if( comparator ) {
            return comparator(k1) <=> comparator(k2)
        }
        else {
            return (k1 as Comparable) <=> (k2 as Comparable)
        }
    }

    protected void visit(int current, Closure closure) {
        assert current >= 0
        assert closure != null

        BTree<K> node = index[current]
        while( node ) {

            if( node.left != -1 )
                visit(node.left, closure)

            closure.call(node.key)

            if( node.right != -1 )
                visit(node.right , closure)
        }
    }

    protected Match<K> nearest(K key) {
        int current = 0
        BTree<K> node = index[current]
        while( node ) {
            int delta = compare(key, node.key)

            if( delta == 0 )
                return new Match<K>(delta:0, node:node, index:current)

            if( delta < 0 && node.left!=-1 ) {
                current = node.left
                node = index[current]
            }
            else if( delta > 0 && node.right!=-1 ) {
                current = node.right
                node = index[current]
            }
            else {
                return new Match<K>(delta:delta, node:node, index:current)
            }

        }
        return null
    }

    protected BTree<K> find(K key) {
        final result = nearest(key)
        result?.delta == 0 ? result.node : null
    }

    @Override
    int size() { store.size() }

    @Override
    boolean isEmpty() { store.isEmpty() }

    @Override
    boolean containsKey(Object key) {
        return index.containsKey(key)
    }

    @Override
    boolean containsValue(Object value) {
        throw new UnsupportedOperationException()
    }

    @Override
    V get(Object key) {
        def node = find((K)key)
        return node ? store[node.pos] : null
    }

    @Override
    V put(K key, V value) {
        assert key

        def result = nearest(key)
        if( result && result.delta == 0 ) {
            final node = result.node
            final old = store[node.pos]
            store[node.pos] = value
            return old
        }

        // save the value
        store[count] = value
        // save a tree node for this value
        index[count] = new BTree<K>( key: key, pos: count )

        if( result ) {
            // update the parent node
            final parent = result.node
            final p = result.index
            if( result.delta < 0 )
                parent.left = count
            else
                parent.right = count
            index[p] = parent
        }

        // increment the count
        count++
        return null
    }

    @Override
    V remove(Object key) {
        throw new UnsupportedOperationException()
    }

    @Override
    void putAll(Map<? extends K, ? extends V> map) {
        for( K key : map.keySet()) {
            put(key, map.get(key))
        }
    }

    @Override
    void clear() {
        store.clear()
        index.clear()
        count=0
    }

    @Override
    Set<K> keySet() {
        return null
    }

    @Override
    Collection<V> values() {
        return null
    }

    @Override
    Set<Map.Entry<K, V>> entrySet() {
        return null
    }

    @Override
    Iterator<K> iterator() { new KeysIterator<K>(this) }
}

@CompileStatic
@ToString(includeNames = true, includePackage = false)
class BTree<T> implements Serializable {
    /*
     * the 'key' of this node
     */
    T key

    /*
     * position where the value is stored
     */
    int pos

    /*
     * position of the left child containing a key less than the one in this node
     */
    int left = -1

    /*
     * position of the right child containing a key greater than the one in this node
     */
    int right = -1
}

@CompileStatic
@ToString(includeNames = true, includePackage = false)
class Match<T> {
    BTree<T> node
    int delta
    int index
}

@CompileStatic
class KeysIterator<T> implements Iterator<T> {

    private SortFileCollector<T,?> collector

    private LinkedList<Integer> stack = new LinkedList<>()

    KeysIterator(SortFileCollector<T,?> owner) {
        this.collector = owner
        deepFirstSearch(0)
    }

    private void deepFirstSearch(int index) {
        def node = getNode(index)
        while( node ) {
            stack.push(index)
            index = node.left
            node = getNode(index)
        }
    }

    private BTree getNode(int index) {
        collector.getIndex().get(index)
    }

    @Override
    boolean hasNext() {
        !stack.isEmpty()
    }

    @Override
    T next() {
        def result = null
        try {
            def index = stack.pop()
            def node = getNode(index)
            if( node ) {
                result = node.key
                deepFirstSearch(node.right)
            }
        }
        catch( NoSuchElementException e ) {
            // ignore and return null
        }

        return result
    }

    @Override
    void remove() {
        throw new UnsupportedOperationException()
    }
}