package nextflow.extension

import nextflow.Session
import org.pf4j.ExtensionPoint
/**
 * Declare the interface for the implementation of a chanel factory
 * extension e.g. Channel.foo.someMethod
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
     * @return Declare the extension factory scope
     */
    String getScope()

}
