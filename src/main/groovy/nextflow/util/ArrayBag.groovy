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
import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.KryoSerializable
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.CompileStatic
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * A bag implementation based on an array list.
 * <p>
 * Note: the main goal of this task is to provide a container class
 * for which the cache hash key is invariant when the order of the
 * items in the container changes. See {@link CacheHelper#hasher(java.lang.Object, nextflow.util.CacheHelper.HashMode)}
 * <p>
 * However the class equals and hashCode methods are the ones provided by
 * the underlying ArrayList i.e. they depend on the content order.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class ArrayBag<E> implements Bag<E>, List<E>, KryoSerializable {

    @Delegate(interfaces = false)
    List target

    ArrayBag() { target = new ArrayList() }

    ArrayBag( int size ) {
        target = new ArrayList(size)
    }

    ArrayBag( Collection items ) {
        target = items != null ? new ArrayList(items) : []
    }

    ArrayBag( Object ... items ) {
        this(items as List)
    }

    @Override
    String toString() {
        InvokerHelper.inspect(this)
    }

//    E getAt( int index )  {
//        target.get(index)
//    }
//
//    void putAt( int index, E value ) {
//        target.set(index, value)
//    }
//
//    @Override
//    int hashCode() {
//        int h = 0;
//        Iterator<E> i = target.iterator();
//        while (i.hasNext()) {
//            E obj = i.next();
//            if (obj != null)
//                h += obj.hashCode();
//        }
//        return h;
//    }
//
//    @Override
//    boolean equals(Object o) {
//        if ( o.is(this) )
//            return true;
//
//        if (!(o instanceof ArrayBag))
//            return false;
//
//        Collection other = ((ArrayBag)o).target
//        if (other.size() != target.size())
//            return false;
//
//        try {
//            return target.containsAll(other);
//        }
//        catch (ClassCastException unused)   {
//            return false;
//        }
//        catch (NullPointerException unused) {
//            return false;
//        }
//    }

    void read (Kryo kryo, Input input) {
        target = kryo.readObject(input,ArrayList)
    }

    void write (Kryo kryo, Output output) {
        kryo.writeObject(output, target)
    }

}
