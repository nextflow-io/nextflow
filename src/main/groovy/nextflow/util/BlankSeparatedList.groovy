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
import groovy.transform.EqualsAndHashCode
import org.codehaus.groovy.runtime.InvokerHelper

/**
 * A list of staged paths, which renders its content just separating
 * the items by a blank space
 */

@CompileStatic
@EqualsAndHashCode
public class BlankSeparatedList implements KryoSerializable {

    // note: this class must NO implement the List interface, otherwise the toString() method is not invoked
    @Delegate(interfaces = false)
    List target

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

    void read (Kryo kryo, Input input) {
        target = kryo.readObject(input,ArrayList)
    }

    void write (Kryo kryo, Output output) {
        kryo.writeObject(output, target)
    }

    def getAt( int index ) {
        target.getAt(index)
    }

    /*
     * this is needed to implement Groovy List extension methods
     */
    def methodMissing(String name, args) {
        InvokerHelper.invokeMethod(target,name,args)
    }

}
