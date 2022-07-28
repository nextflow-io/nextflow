package nextflow.extension

import groovy.transform.PackageScope
import nextflow.Session
import org.pf4j.ExtensionPoint
/**
 * Define a function factory extension. A plugin can provide a function factory
 * by implementing this interface. Public methods in such class are accessible
 * as function extensions over the declared scope.
 *
 * For example having the `foo` scope and the `bar` method. It will be possible
 * to invoke `bar()`
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class FunctionExtensionPoint implements ExtensionPoint {

    private boolean initialised

    @PackageScope
    synchronized void checkInit(Session session) {
        if( !initialised ) {
            initialised = true
        }
    }

}
