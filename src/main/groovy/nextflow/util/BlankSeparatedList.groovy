package nextflow.util
/**
 * A list of staged paths, which renders its content just separating
 * the items by a blank space
 */

public class BlankSeparatedList {

    // note: due to a bug with the groovy runtime, the class must no implement
    // the List interface, otherwise the toString() method is not invoked (Groovy 2.1.7)
    @Delegate(interfaces = false)
    List target

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

}
