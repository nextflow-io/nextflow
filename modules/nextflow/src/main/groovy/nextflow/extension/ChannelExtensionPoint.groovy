package nextflow.extension

import groovy.transform.PackageScope
import nextflow.Session
import org.pf4j.ExtensionPoint
/**
 * Define a channel factory extension. A plugin can extend a channel factory
 * by implementing this interface. Public methods in such class are accessible
 * as channel extensions over the declared scope.
 *
 * For example having the `foo` scope and the `bar` method. It will be possible
 * to invoke `channel.foo.bar()` 
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class ChannelExtensionPoint implements ExtensionPoint {

    private boolean initialised

    @PackageScope
    synchronized void checkInit(Session session) {
        if( !initialised ) {
            init(session)
            initialised = true
        }
    }


    /**
     * Channel factory initialization. This method is invoked one and only once before
     * the before target extension method is called.
     *
     * @param session The current nextflow {@link Session}
     */
    abstract protected void init(Session session)

}
