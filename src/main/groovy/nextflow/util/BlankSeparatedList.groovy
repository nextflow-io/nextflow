package nextflow.util
import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.KryoSerializable
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.transform.EqualsAndHashCode
/**
 * A list of staged paths, which renders its content just separating
 * the items by a blank space
 */

@EqualsAndHashCode
public class BlankSeparatedList implements KryoSerializable {

    // note: due to a bug with the groovy runtime, the class must NO implement
    // the List interface, otherwise the toString() method is not invoked (Groovy 2.1.7)
    @Delegate(interfaces = false)
    List target

    public BlankSeparatedList() { }

    public BlankSeparatedList( Collection items ) {
        target = items != null ? new ArrayList(items) : []
    }

    public BlankSeparatedList( Object ... items ) {
        this(items as List)
    }

    @Override
    public String toString() {
        target.join(' ')
    }


    public void read (Kryo kryo, Input input) {
        target = kryo.readObject(input,ArrayList)
    }

    public void write (Kryo kryo, Output output) {
        kryo.writeObject(output, target)
    }

}
