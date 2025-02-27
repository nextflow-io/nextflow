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

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.KryoSerializable
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
/**
 * A list of staged paths, which renders its content just separating
 * the items by a blank space
 */

@CompileStatic
@EqualsAndHashCode
class BlankSeparatedList implements KryoSerializable, PathEscapeAware, List {

    private List target

    // note: this constructor is needed by kryo serialization
    private BlankSeparatedList() { }

    BlankSeparatedList( Collection items ) {
        target = items != null ? new ArrayList(items) : []
    }

    BlankSeparatedList( Object ... items ) {
        this(items as List)
    }

    @Override
    String toString() {
        target.join(' ')
    }

    String toStringEscape() {
        def result = new StringBuilder()
        for( int i=0; i<target.size(); i++ ) {
            if( i ) result.append(' ')
            result.append(Escape.path(target.get(i).toString()))
        }
        return result.toString()
    }

    void read (Kryo kryo, Input input) {
        target = kryo.readObject(input,ArrayList)
    }

    void write (Kryo kryo, Output output) {
        kryo.writeObject(output, target)
    }

    def getAt( int index ) {
        target.getAt(index)
    }

    @Override
    int size() {
        return target.size()
    }

    @Override
    boolean isEmpty() {
        return target.isEmpty()
    }

    @Override
    boolean contains(Object o) {
        return target.contains(o)
    }

    @Override
    Iterator iterator() {
        return target.iterator()
    }

    @Override
    Object[] toArray() {
        return target.toArray()
    }

    @Override
    boolean add(Object object) {
        return target.add(object)
    }

    @Override
    boolean remove(Object o) {
        return target.remove(o)
    }

    @Override
    boolean addAll(Collection c) {
        return target.addAll(c)
    }

    @Override
    boolean addAll(int index, Collection c) {
        return target.addAll(index, c)
    }

    @Override
    void clear() {
        target.clear()
    }

    @Override
    Object get(int index) {
        return target.get(index)
    }

    @Override
    Object set(int index, Object element) {
        return target.set(index, element)
    }

    @Override
    void add(int index, Object element) {
        target.add(index, element)
    }

    @Override
    Object remove(int index) {
        return target.remove(index)
    }

    @Override
    int indexOf(Object o) {
        return target.indexOf(o)
    }

    @Override
    int lastIndexOf(Object o) {
        return target.lastIndexOf(o)
    }

    @Override
    ListIterator listIterator() {
        return target.listIterator()
    }

    @Override
    ListIterator listIterator(int index) {
        return target.listIterator(index)
    }

    @Override
    List subList(int fromIndex, int toIndex) {
        return target.subList(fromIndex, toIndex)
    }

    @Override
    boolean retainAll(Collection c) {
        return target.retainAll(c)
    }

    @Override
    boolean removeAll(Collection c) {
        return target.remove(c)
    }

    @Override
    boolean containsAll(Collection c) {
        return target.containsAll(c)
    }

    @Override
    Object[] toArray(Object[] a) {
        return target.toArray(a)
    }
}
