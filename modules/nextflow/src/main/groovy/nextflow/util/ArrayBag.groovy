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
