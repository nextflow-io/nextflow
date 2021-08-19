package nextflow.extension

import java.lang.reflect.Method

import groovy.runtime.metaclass.ChannelFactory
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.plugin.Plugins
import nextflow.plugin.Scoped
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Object holding a set of {@link ChannelExtensionPoint} instances
 * for a given scope.
 *
 * It collects all channel factory extension classes provided by plugins
 * for a scope context.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Canonical
class ChannelFactoryImpl implements ChannelFactory {

    @PackageScope
    static List<ChannelExtensionPoint> allExtensions

    /**
     * The scope of implemented by this channel factory holder. For example
     * having extension `Channel.foo.something()` the scope is `foo` 
     */
    private String scope

    /**
     * The set of channel extensions held by this scope
     */
    private ChannelExtensionPoint target

    ChannelFactoryImpl(String scope, ChannelExtensionPoint extensionClass) {
        this.scope = scope
        this.target = extensionClass
    }


    /**
     * Implements extension method invocation login
     *
     * @param methodName The name of the extension method to be invoked
     * @param args The array holding the arguments method arguments
     * @return The extension method result
     */
    private Object invoke0(String methodName, Object[] args) {
        final meta = target.metaClass.getMetaMethod(methodName, args)
        if( meta && meta.isPublic() ) {
            target.checkInit((Session)Global.session)
            final Method method = target.getClass().getMethod(methodName, meta.getNativeParameterTypes())
            // or -- owner.metaClass.invokeMethod(target, methodName, args)
            return method.invoke(target, args)
        }
        throw new MissingFactoryMethodException(this.scope, methodName, args)
    }

    /**
     * Invokes the extension method with the given arguments
     *
     * @param methodName The name of the extension method to be invoked
     * @param args The array holding the arguments method arguments
     * @return The extension method result
     */
    @Override
    Object invokeExtensionMethod(String method, Object[] args) {
        def result = invoke0(method,args)
        NodeMarker.addSourceNode("channel.${this.scope}.${method}", result)
        return result
    }

    /**
     * Create a {@link ChannelFactoryImpl} object for the given scope
     *
     * @param scope
     * @return The {@link ChannelFactoryImpl} instance holding the extensions matching the specified scope
     */
    static ChannelFactoryImpl create(String scope) {
        final all = findAllExtensions()
        log.debug "Loading channel factory extensions: $all"
        for( ChannelExtensionPoint it : all ) {
            final annot = it.getClass().getAnnotation(Scoped)
            if( annot && annot.value()==scope )
                return new ChannelFactoryImpl(scope, it)
        }
        return null
    }

    static protected List<ChannelExtensionPoint> findAllExtensions() {
        if( allExtensions==null )
            allExtensions = Plugins.getScopedExtensions(ChannelExtensionPoint).toList()
        return allExtensions
    }

    /**
     * Customized missing method  extension
     */
    static class MissingFactoryMethodException extends MissingMethodException {

        private String scope

        MissingFactoryMethodException(String scope, String method, Object[] args) {
            super(method, Channel, args)
            this.scope = scope
        }

        @Override
        String getMessage() {
            return "No signature of method: " +
                    "Channel.${this.scope}.${method}" +
                    "() is applicable for argument types: (" +
                    InvokerHelper.toTypeString(arguments, 60) +
                    ") values: " +
                    InvokerHelper.toArrayString(arguments, 60, true)
        }

        String toString() {
            return getMessage()
        }
    }
}
