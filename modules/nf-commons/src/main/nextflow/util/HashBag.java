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

package nextflow.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

import nextflow.script.types.Bag;

/**
 * Implementation of a bag based on a frequency map.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class HashBag<E> implements Bag<E> {

    private Map<E, Integer> counts;

    public HashBag() {
        this(10);
    }

    public HashBag(int initialCapacity) {
        counts = new HashMap<>(initialCapacity);
    }

    public HashBag(Collection<? extends E> c) {
        this(c.size());
        addAll(c);
    }

    /**
     * Ensures that this collection contains the specified element.
     *
     * @param e
     */
    @Override
    public boolean add(E e) {
        var v = counts.computeIfAbsent(e, (k) -> 0);
        counts.put(e, v + 1);
        return true;
    }

    /**
     * Adds all of the elements in the specified collection to this collection.
     *
     * @param c
     */
    @Override
    public boolean addAll(Collection<? extends E> c) {
        for( var e : c )
            add(e);
        return true;
    }

    /**
     * Removes all of the elements from this collection.
     */
    @Override
    public void clear() {
        counts.clear();
    }

    /**
     * Returns true if this collection contains the specified element.
     *
     * @param o
     */
    @Override
    public boolean contains(Object o) {
        return counts.containsKey(o);
    }

    /**
     * Returns true if this collection contains all of the elements in the specified collection.
     *
     * @param c
     */
    @Override
    public boolean containsAll(Collection<?> c) {
        for( var e : c ) {
            if( !contains(e) )
                return false;
        }
        return true;
    }

    /**
     * Compares the specified object with this collection for equality.
     *
     * @param o
     */
    @Override
    public boolean equals(Object o) {
        return o instanceof HashBag that && this.counts.equals(that.counts);
    }

    /**
     * Returns the hash code value for this collection.
     */
    @Override
    public int hashCode() {
        return counts.hashCode();
    }

    /**
     * Returns true if this collection contains no elements.
     */
    @Override
    public boolean isEmpty() {
        return counts.isEmpty();
    }

    /**
     * Returns an iterator over the elements in this collection.
     */
    @Override
    public Iterator<E> iterator() {
        return new Iterator<E>() {
            private Iterator<Map.Entry<E,Integer>> itor = counts.entrySet().iterator();
            private E currentElement;
            private int currentCount;

            @Override
            public boolean hasNext() {
                return currentCount > 0 || itor.hasNext();
            }

            @Override
            public E next() {
                if( currentCount == 0 ) {
                    var entry = itor.next();
                    currentElement = entry.getKey();
                    currentCount = entry.getValue();
                }
                currentCount--;
                return currentElement;
            }
        };
    }

    /**
     * Removes a single instance of the specified element from this collection, if it is present.
     *
     * @param o
     */
    @Override
    public boolean remove(Object o) {
        var count = counts.get(o);
        if( count == null )
            return false;
        if( count.intValue() == 1 )
            counts.remove(o);
        else
            counts.put((E) o, count.intValue() - 1);
        return true;
    }

    /**
     * Removes all of this collection's elements that are also contained in the specified collection.
     *
     * @param c
     */
    @Override
    public boolean removeAll(Collection<?> c) {
        var changed = false;
        for( var e : c ) {
            if( remove(e) )
                changed = true;
        }
        return changed;
    }

    /**
     * Retains only the elements in this collection that are contained in the specified collection.
     *
     * @param c
     */
    @Override
    public boolean retainAll(Collection<?> c) {
        var changed = false;
        for( var e : new HashSet<>(counts.keySet()) ) {
            if( !c.contains(e) ) {
                counts.remove(e);
                changed = true;
            }
        }
        return changed;
    }

    /**
     * Returns the number of elements in this collection.
     */
    @Override
    public int size() {
        var result = 0;
        for( var count : counts.values() )
            result += count;
        return result;
    }

    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
        var builder = new StringBuilder();
        builder.append('[');
        var itor = iterator();
        while( itor.hasNext() ) {
            builder.append(itor.next().toString());
            if( itor.hasNext() )
                builder.append(", ");
        }
        builder.append(']');
        return builder.toString();
    }

}
