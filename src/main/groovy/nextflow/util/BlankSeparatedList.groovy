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
public class BlankSeparatedList implements KryoSerializable {

    // note: due to a bug with the groovy runtime, the class must NO implement
    // the List interface, otherwise the toString() method is not invoked (Groovy 2.1.7)
    @Delegate(interfaces = false)
    List target

    BlankSeparatedList() { }

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

}
