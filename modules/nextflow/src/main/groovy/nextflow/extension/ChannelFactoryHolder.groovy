package nextflow.extension

import java.lang.reflect.Method
import java.util.concurrent.ConcurrentHashMap

import groovy.runtime.metaclass.ChannelFactoryExtension
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker
import nextflow.plugin.Plugins
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Object holding a set of {@link ChannelExtensionPoint} instances
 * for a given scope. In practical term it collect all channel factory extensions
 * for a scope.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Canonical
class ChannelFactoryHolder implements ChannelFactoryExtension {

    /**
     * The scope of implemented by this channel factory holder. For example
     * having extension `Channel.foo.something()` the scope is `foo` 
     */
    private String scope

    /**
     * The set of channel extensions held by this scope
     */
    private Set<ChannelExtensionPoint> exts = new HashSet<>()

    private ConcurrentHashMap<ChannelExtensionPoint,Boolean> markInitialized = new ConcurrentHashMap<>()

    ChannelFactoryHolder(String scope, Collection<ChannelExtensionPoint> extensions) {
        this.scope = scope
        this.exts = new HashSet<>(extensions)
    }

    /*
     * Invoke the init method exactly once
     */
    protected void checkInit(ChannelExtensionPoint target) {
        def wasInit = markInitialized.putIfAbsent(target, true)
        if( !wasInit )
            target.init(Global.session as Session)
    }

    private Object invoke0(String methodName, Object[] args) {
        for( int i=0; i<exts.size(); i++ ) {
            final target = exts[i]
            final meta = target.metaClass.getMetaMethod(methodName, args)
            if( meta ) {
                checkInit(target)
                final Method method = target.getClass().getMethod(methodName, meta.getNativeParameterTypes())
                // or -- owner.metaClass.invokeMethod(owner, methodName, args)
                return method.invoke(target, args)
            }
        }

        throw new MissingFactoryMethodException(this.scope, methodName, args)
    }
    
    @Override
    Object invokeExtensionMethod(String method, Object[] args) {
        def result = invoke0(method,args)
        NodeMarker.addSourceNode("Channel.${this.scope}.${method}", result)
        return result
    }

    /**
     * Create a {@link ChannelFactoryHolder} object for the given scope
     *
     * @param scope
     * @return The {@link ChannelFactoryHolder} instance holding the extensions matching the specified scope
     */
    static ChannelFactoryHolder create(String scope) {
        final all = Plugins.getExtensions(ChannelExtensionPoint)
        log.debug "Loading factory extensions: $all"
        final matchingType = []
        for( ChannelExtensionPoint ext : all ) {
            if( ext.scope == scope ) {
                matchingType.add(ext)
            }
        }
        return new ChannelFactoryHolder(scope, matchingType)
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
