package nextflow.extension

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
interface ChannelExtensionPoint extends ExtensionPoint {

    /**
     * Channel factory initialization. This method is invoked one and only once before
     * the before target extension method is called.
     *
     * @param session The current nextflow {@link Session}
     */
    void init(Session session)

    /**
     * @return Declare the channel factory scope
     */
    String getScope()

}
