/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.PoisonPill
import groovyx.gpars.dataflow.stream.DataflowStreamReadAdapter

/**
 * This class wrap a channel behaving like a {@code List}
 */
class ChannelListProxy implements List<Object> {

    @Lazy
    private List target = { getRealList() }()

    final private DataflowReadChannel channel

    ChannelListProxy( DataflowReadChannel channel ) {
        this.channel = channel
    }

    private List getRealList() {

        List result = []
        switch(channel) {
            case DataflowVariable:
                def value = channel.getVal()
                result.add( value )
                break

            case DataflowQueue:
            case DataflowStreamReadAdapter:
                while(true) {
                    def value = channel.getVal()
                    if( value instanceof PoisonPill) {
                        break
                    }
                    result.add(value)
                }
                break

            default:
                throw new IllegalArgumentException("Invalid channel type: ${channel?.class?.name}")
        }

        return result
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
    Iterator<Object> iterator() {
        return target.iterator()
    }

    @Override
    Object[] toArray() {
        return target.toArray()
    }

    @Override
    def <T> T[] toArray(T[] a) {
        return target.toArray(a)
    }

    @Override
    boolean add(Object e) {
        return target.add(e)
    }

    @Override
    boolean remove(Object o) {
        return target.remove(o)
    }

    @Override
    Object remove(int index) {
        return target.remove(index)
    }

    @Override
    boolean containsAll(Collection<?> c) {
        return target.containsAll(c)
    }

    @Override
    boolean addAll(Collection<?> c) {
        return target.addAll(c)
    }

    @Override
    boolean addAll(int index, Collection<?> c) {
        return target.addAll(index, c)
    }

    @Override
    boolean removeAll(Collection<?> c) {
        return target.removeAll(c)
    }

    @Override
    boolean retainAll(Collection<?> c) {
        return target.retainAll(c)
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
    int indexOf(Object o) {
        return target.indexOf(o)
    }

    @Override
    int lastIndexOf(Object o) {
        return target.lastIndexOf(o)
    }

    @Override
    ListIterator<Object> listIterator() {
        return target.listIterator()
    }

    @Override
    ListIterator<Object> listIterator(int index) {
        return  target.listIterator(index)
    }

    @Override
    List<Object> subList(int fromIndex, int toIndex) {
        return target.subList(fromIndex, toIndex)
    }
}