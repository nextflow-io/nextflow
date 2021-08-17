package nextflow.extension

import java.lang.reflect.Method
import java.util.concurrent.ConcurrentHashMap

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
    private Set<ChannelExtensionPoint> extensions = new HashSet<>()

    private ConcurrentHashMap<ChannelExtensionPoint,Boolean> markInitialized = new ConcurrentHashMap<>()

    ChannelFactoryImpl(String scope, Collection<ChannelExtensionPoint> extensions) {
        this.scope = scope
        this.extensions = new HashSet<>(extensions)
    }

    /*
     * Invoke the init method exactly once
     */
    protected void checkInit(ChannelExtensionPoint target) {
        final wasInit = markInitialized.putIfAbsent(target, true)
        if( !wasInit )
            target.init(Global.session as Session)
    }

    private Object invoke0(String methodName, Object[] args) {
        for(int i=0; i<extensions.size(); i++ ) {
            final target = extensions[i]
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
     * Create a {@link ChannelFactoryImpl} object for the given scope
     *
     * @param scope
     * @return The {@link ChannelFactoryImpl} instance holding the extensions matching the specified scope
     */
    static ChannelFactoryImpl create(String scope) {
        final all = findAllExtensions()
        log.debug "Loading channel factory extensions: $all"
        final matchingClasses = new ArrayList(10)
        for( ChannelExtensionPoint it : all ) {
            if( it.getScope() == scope ) {
                matchingClasses.add(it)
            }
        }
        return matchingClasses ? new ChannelFactoryImpl(scope, matchingClasses) : null
    }

    static protected List<ChannelExtensionPoint> findAllExtensions() {
        if( allExtensions==null )
            allExtensions = Plugins.getExtensions(ChannelExtensionPoint)
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
