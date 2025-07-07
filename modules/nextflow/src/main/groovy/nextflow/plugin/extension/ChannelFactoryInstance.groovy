/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.plugin.extension

import groovy.runtime.metaclass.ChannelFactory
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Global
import nextflow.Session
import nextflow.dag.NodeMarker

/**
 * Object holding a set of {@link PluginExtensionPoint} instances
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
class ChannelFactoryInstance implements ChannelFactory {

    /**
     * The set of channel extensions held by this scope
     */
    private PluginExtensionPoint target

    ChannelFactoryInstance(PluginExtensionPoint extensionClass) {
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
            final method = target.getClass().getMethod(methodName, meta.getNativeParameterTypes())
            // fix issue casting issue when argument is a GString but target method expects a String
            for( int i=0; i<args.length; i++ )
                if( args[i] instanceof GString && method.getParameterTypes()[i].isAssignableFrom(String) ) args[i] = args[i].toString()
            // finally invoke it
            return method.invoke(target, args)
        }
        throw new MissingMethodException(methodName, Channel, args)
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
        NodeMarker.addSourceNode("channel.${method}", result)
        return result
    }

}
