/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.extension

import groovy.transform.MapConstructor

import java.lang.reflect.Modifier

import groovy.runtime.metaclass.ExtensionProvider
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.Session
import nextflow.plugin.Plugins
import nextflow.script.ChannelOut
/**
 * Manage channel extensions and dispatch method invocations
 * to target class implementing the extension logic
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ChannelExtensionProvider implements ExtensionProvider {

    private static ChannelExtensionProvider instance

    private Session getSession() { Global.getSession() as Session }

    /**
     * Hold all available operator extensions. The key represent the name of operator how it's expected
     * to be invoked (it can be the alias name). The value is an object holding the real method name
     * and the target object on which it operator will be invoked
     */
    final private Map<String,PluginExtensionMethod> operatorExtensions = new HashMap<>()

    /**
     * Hold all available factory extensions. The key represent the name of operator how it's expected
     * to be invoked (it can be the alias name). The value is an object holding the real method name
     * and the target object on which it operator will be invoked
     */
    final private Map<String,PluginExtensionMethod> factoryExtensions = new HashMap<>()

    private List<ChannelExtensionPoint> channelExtensionPoints

    private Set<String> OPERATOR_NAMES

    static ChannelExtensionProvider INSTANCE() {
        if( instance != null )
            return instance
        return instance = new ChannelExtensionProvider().install()
    }

    static void reset() {
        instance = null
    }

    ChannelExtensionProvider install() {
        // add default operators
        final defaultOps = loadDefaultOperators()
        log.trace "Dataflow default extension methods: ${defaultOps.sort().join(',')}"
        OPERATOR_NAMES = Collections.unmodifiableSet(defaultOps)
        // configure as global instance
        return instance = this
    }

    /**
     * Load all operators defined by nextflow
     * @return The set of operator names
     */
    private Set<String> loadDefaultOperators() {
        final result = getDeclaredExtensionMethods0(OperatorEx.class)
        for( String it : result )
            operatorExtensions.put(it, new PluginExtensionMethod(method: it, target: OperatorEx.instance))
        return result
    }

    /**
     * Load all extension method declared by the specified plugin Id
     *
     * @param pluginId The Id of the plugin from where the extension methods need to be loaded
     * @param includedNames The map of extension method as provided in the `include` declaration.
     *      The key represent the real method name and the value the name alias the method needs to
     *      be referenced in the script. If the alias is not provided the key == value.
     * @return
     *      The class itself to allow method chaining
     */
    ChannelExtensionProvider loadPluginExtensionMethods(String pluginId, Map<String, String> includedNames){
        final ext= findPluginExtensionMethods(pluginId)
        if( ext ) {
            loadPluginExtensionMethods(ext, includedNames)
        }
        return instance = this
    }

    protected ChannelExtensionProvider loadPluginExtensionMethods(ChannelExtensionPoint ext, Map<String, String> includedNames){
        // find all operators defined in the plugin
        final definedOperators= getDeclaredExtensionMethods0(ext.getClass())
        // final all factories defined in the plugin
        final definedFactories= getDeclaredFactoryExtensionMethods0(ext.getClass())
        for(Map.Entry<String,String> entry : includedNames ) {
            String realName = entry.key
            String aliasName = entry.value
            final reference = operatorExtensions.get(aliasName)
            if( reference ){
                throw new IllegalStateException("Operator '$aliasName' conflict - it's defined by plugin ${reference.target.class.name}")
            }
            Object existing = operatorExtensions.get(aliasName)
            if (existing.is(OperatorEx.instance)) {
                throw new IllegalStateException("Operator '$realName' is already defined as a built-in operator - Offending plugin class: $ext")
            }
            else if (existing != null) {
                if( existing.getClass().getName() != ext.getClass().getName() ) {
                    throw new IllegalStateException("Operator '$realName' conflict - it's defined by plugin ${existing.getClass().getName()} and ${ext.getClass().getName()}")
                }
            }
            if( definedOperators.contains(realName) ) {
                OPERATOR_NAMES = Collections.unmodifiableSet(OPERATOR_NAMES + [aliasName])
                operatorExtensions.put(aliasName, new PluginExtensionMethod(method:realName, target:ext))
            }
            else if( definedFactories.contains(realName) ){
                ChannelFactoryInstance factoryInstance = new ChannelFactoryInstance(ext)
                factoryExtensions.put(aliasName, new PluginExtensionMethod(method:realName, target:factoryInstance))
            }
            else{
                throw new IllegalStateException("Operator '$realName' it isn't defined by plugin ${existing.getClass().getName()}")
            }
        }
        return instance = this
    }

    static private Set<String> getDeclaredExtensionMethods0(Class clazz) {
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) ) continue
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) ) continue
            // operator extension method must have a dataflow read channel type as first argument
            def params=handle.getParameterTypes()
            if( params.length>0 && isReadChannel(params[0]) )
                result.add(handle.name)
        }
        return result
    }

    static private Set<String> getDeclaredFactoryExtensionMethods0(Class clazz) {
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) ) continue
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) ) continue
            // factory extension method must have a dataflow write channel type as return
            def returnType =handle.getReturnType()
            if( isWriteChannel(returnType) )
                result.add(handle.name)
        }
        return result
    }

    static boolean isReadChannel(Class clazz) {
        DataflowReadChannel.class.isAssignableFrom(clazz)
    }

    static boolean isWriteChannel(Class clazz) {
        DataflowWriteChannel.class.isAssignableFrom(clazz)
    }

    Set<String> operatorNames() { OPERATOR_NAMES }

    boolean isExtensionMethod(Object obj, String name) {
        if( obj instanceof DataflowReadChannel || obj instanceof DataflowBroadcast || obj instanceof ChannelOut ) {
            return OPERATOR_NAMES.contains(name)
        }
        return false
    }

    Object invokeExtensionMethod(Object channel, String method, Object[] args) {
        final target = operatorExtensions.get(method)
        if( target==null )
            throw new IllegalStateException("Missing target class for operator '$method'")
        method = operatorExtensions.get(method)?.method
        if( target.target instanceof ChannelExtensionPoint )
            ((ChannelExtensionPoint)target.target).checkInit(getSession())
        new OpCall(target.target,channel,method,args).call()
    }

    def invokeFactoryExtensionMethod(String name, Object[] args){
        if( factoryExtensions.containsKey(name) ){
            def reference = factoryExtensions.get(name)
            def factory = (ChannelFactoryInstance)reference.target
            return factory.invokeExtensionMethod(reference.method, args)
        }
    }

    protected ChannelExtensionPoint findPluginExtensionMethods(String pluginId) {
        Plugins.getExtensionsInPluginId(ChannelExtensionPoint, pluginId)?.first()
    }

    @Deprecated
    static void reloadExtensionPoints() {
        if( !instance )
            return
        instance.channelExtensionPoints=null
        instance.operatorExtensions.clear()
        instance.install()
    }

    /**
     * Hold a reference to a extension method provided by a Nextflow plugin
     */
    @MapConstructor
    class PluginExtensionMethod {
        /**
         * The name of the method that needs to be invoked
         */
        String method

        /**
         * The target object on which the method is going to be invoked
         */
        Object target
    }
}
