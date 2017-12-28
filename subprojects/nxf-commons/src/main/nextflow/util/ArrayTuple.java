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

package nextflow.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 * Provides a basic tuple implementation extending an {@link ArrayList}
 * and not allowing any content change operation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class ArrayTuple<E> extends ArrayList<E> {

    private static final long serialVersionUID = - 4765828600345948947L;

    public ArrayTuple() {}

    public ArrayTuple(Collection<E> other) {
        super(other);
    }

    public ArrayTuple(int initialCapacity) {
        super(initialCapacity);
    }

    public boolean add(E e) {
        throw new UnsupportedOperationException();
    }

    public boolean remove(Object o) {
        throw new UnsupportedOperationException();
    }

    public boolean addAll(Collection<? extends E> coll) {
        throw new UnsupportedOperationException();
    }

    public boolean removeAll(Collection<?> coll) {
        throw new UnsupportedOperationException();
    }

    public boolean retainAll(Collection<?> coll) {
        throw new UnsupportedOperationException();
    }

    public void clear() {
        throw new UnsupportedOperationException();
    }

    public E set(int index, E element) {
        throw new UnsupportedOperationException();
    }

    public void add(int index, E element) {
        throw new UnsupportedOperationException();
    }

    public E remove(int index) {
        throw new UnsupportedOperationException();
    }

    public boolean addAll(int index, Collection<? extends E> c) {
        throw new UnsupportedOperationException();
    }

    public Iterator<E> iterator() {
        return new Iterator<E>() {
            private final Iterator<? extends E> i = ArrayTuple.super.iterator();

            public boolean hasNext() {return i.hasNext(); }
            public E next()          {return i.next(); }
            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public ListIterator<E> listIterator(final int index) {
        return new ListIterator<E>() {
            private final ListIterator<? extends E> i = ArrayTuple.super.listIterator(index);

            public boolean hasNext()     {return i.hasNext();}
            public E next()              {return i.next();}
            public boolean hasPrevious() {return i.hasPrevious();}
            public E previous()          {return i.previous();}
            public int nextIndex()       {return i.nextIndex();}
            public int previousIndex()   {return i.previousIndex();}

            public void remove() {
                throw new UnsupportedOperationException();
            }
            public void set(E e) {
                throw new UnsupportedOperationException();
            }
            public void add(E e) {
                throw new UnsupportedOperationException();
            }
        };
    }

    public List<E> subList(int fromIndex, int toIndex) {
        return new ArrayTuple<E>(super.subList(fromIndex, toIndex));
    }


}