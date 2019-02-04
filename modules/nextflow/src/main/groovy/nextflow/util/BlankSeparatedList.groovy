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
