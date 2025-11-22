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

import java.lang.reflect.Modifier

import groovy.runtime.metaclass.ExtensionProvider
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.extension.OpCall
import nextflow.extension.OperatorImpl
import nextflow.plugin.Plugins
import nextflow.script.ChannelOut
import nextflow.script.FunctionDef
import nextflow.script.ScriptMeta
import nextflow.util.TestOnly

/**
 * Manage channel extensions and dispatch method invocations
 * to target class implementing the extension logic
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class PluginExtensionProvider implements ExtensionProvider {

    private static PluginExtensionProvider instance

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

    private Set<String> OPERATOR_NAMES

    static PluginExtensionProvider INSTANCE() {
        if( instance != null )
            return instance
        return instance = new PluginExtensionProvider().install()
    }

    @TestOnly
    static void reset() {
        instance = null
    }

    PluginExtensionProvider install() {
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
        final result = getDeclaredOperatorExtensionMethods0(OperatorImpl.class, true)
        for( String it : result )
            operatorExtensions.put(it, new PluginExtensionMethod(method: it, target: OperatorImpl.instance))
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
    PluginExtensionProvider loadPluginExtensionMethods(String pluginId, Map<String, String> includedNames){
        final extensions= Plugins.getExtensionsInPluginId(PluginExtensionPoint, pluginId)
        if( !extensions )
            throw new AbortOperationException("Plugin '$pluginId' does not implement any extension point")
        if( extensions.size()>1 )
            throw new AbortOperationException("Plugin '$pluginId' implements more than one extension point: ${extensions.collect(it -> it.class.getSimpleName()).join(',')}")
        loadPluginExtensionMethods(pluginId,extensions.first(), includedNames)
        return instance = this
    }

    protected PluginExtensionProvider loadPluginExtensionMethods(String pluginId, PluginExtensionPoint ext, Map<String, String> includedNames){
        // find all operators defined in the plugin
        final definedOperators= getDeclaredOperatorExtensionMethods0(ext.getClass())
        // find all factories defined in the plugin
        final definedFactories= getDeclaredFactoryExtensionMethods0(ext.getClass())
        // find all functions defined in the plugin
        final definedFunctions= getDeclaredFunctionsExtensionMethods0(ext.getClass())
        for( Map.Entry<String,String> entry : includedNames ) {
            String realName = entry.key
            String aliasName = entry.value
            // check if it has already been included
            final existing = operatorExtensions.get(aliasName)
            if (existing.is(OperatorImpl.instance)) {
                throw new IllegalStateException("Operator '$realName' is already defined as a built-in operator - Offending plugin '$pluginId'")
            }
            else if (existing != null) {
                if( existing.getClass().getName() != ext.getClass().getName() ) {
                    throw new IllegalStateException("Operator '$realName' conflict - it's defined by plugin ${pluginId} and ${existing.pluginId}")
                }
            }
            if( definedOperators.contains(realName) ) {
                OPERATOR_NAMES = Collections.unmodifiableSet(OPERATOR_NAMES + [aliasName])
                operatorExtensions.put(aliasName, new PluginExtensionMethod(method:realName, target:ext, pluginId:pluginId))
            }
            else if( definedFactories.contains(realName) ){
                ChannelFactoryInstance factoryInstance = new ChannelFactoryInstance(ext)
                factoryExtensions.put(aliasName, new PluginExtensionMethod(method:realName, target:factoryInstance, pluginId:pluginId))
            }
            else if( definedFunctions.contains(realName) ){
                FunctionDef functionDef = new FunctionDef(ext, realName, aliasName )
                meta.addDefinition(functionDef)
            }
            else{
                throw new IllegalStateException("Extension '$realName' it isn't defined by plugin ${pluginId}")
            }
        }
        // initialise the plugin session
        ext.checkInit((Session)Global.session)
        return instance = this
    }

    static private Set<String> getDeclaredOperatorExtensionMethods0(Class clazz, boolean internal=false) {
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            if( result.contains(handle.name))
                continue
            // in a future only annotated methods will be imported
            if( !internal && handle.isAnnotationPresent(Operator)) {
                if( !Modifier.isPublic(handle.getModifiers()) )
                    throw new IllegalStateException("Operator extension '$handle.name' in `$clazz.name` should be declared public")
                if( Modifier.isStatic(handle.getModifiers()) )
                    throw new IllegalStateException("Operator extension '$handle.name' in `$clazz.name` cannot be not declared as a static method")
                final params=handle.getParameterTypes()
                if( params.length == 0 || !isReadChannel(params[0]) ) {
                    throw new IllegalStateException("Operator extension '$handle.name' in `$clazz.name` has not a valid signature")
                }
                result.add(handle.name)
                continue
            }

            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) )
                continue
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) )
                continue
            // operator extension method must have a dataflow read channel type as first argument
            final params=handle.getParameterTypes()
            if( params.length>0 && isReadChannel(params[0]) ) {
                if( !internal )
                    log.warn("Operator extension `$handle.name` in `$clazz.name` should be marked with the '@Operator' annotation")
                result.add(handle.name)
            }
        }
        return result
    }

    static private Set<String> getDeclaredFactoryExtensionMethods0(Class clazz) {
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            // skip duplicates
            if( result.contains(handle.name)) continue
            // in a future only annotated methodS will be imported
            if( handle.isAnnotationPresent(Factory)) {
                if( !Modifier.isPublic(handle.getModifiers()) )
                    throw new IllegalStateException("Factory extension '$handle.name' in `$clazz.name` should be declared public")
                if( Modifier.isStatic(handle.getModifiers()) )
                    throw new IllegalStateException("Factory extension '$handle.name' in `$clazz.name` cannot be not declared as a static method")
                final returnType = handle.getReturnType()
                if( !isWriteChannel(returnType) )
                    throw new IllegalStateException("Factory extension '$handle.name' in `$clazz.name` has not a valid signature")
                result.add(handle.name)
                continue
            }
            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) ) continue
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) ) continue
            // factory extension method must have a dataflow write channel type as return
            final params=handle.getParameterTypes()
            final returnType = handle.getReturnType()
            if( isWriteChannel(returnType) && (!params || !isReadChannel(params[0])) ) {
                log.warn("Factory extension '$handle.name' in `$clazz.name` should be marked with the '@Factory' annotation")
                result.add(handle.name)
            }
        }
        return result
    }

    static private Set<String>getDeclaredFunctionsExtensionMethods0(Class clazz){
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            // skip duplicates
            if( result.contains(handle.name))
                continue
            // custom functions must to be annotated with @Function
            if( !handle.isAnnotationPresent(Function))
                continue
            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) )
                throw new IllegalStateException("Function extension '$handle.name' in `$clazz.name` should be declared public")
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) )
                throw new IllegalStateException("Function extension '$handle.name' in `$clazz.name` cannot be not declared as a static method")
            result.add(handle.name)
        }
        return result
    }

    @PackageScope
    ScriptMeta getMeta() { ScriptMeta.current() }

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
        if( target.target instanceof PluginExtensionPoint )
            ((PluginExtensionPoint)target.target).checkInit(getSession())
        new OpCall(target.target,channel,method,args).call()
    }

    def invokeFactoryExtensionMethod(String name, Object[] args){
        if( factoryExtensions.containsKey(name) ){
            def reference = factoryExtensions.get(name)
            def factory = (ChannelFactoryInstance)reference.target
            return factory.invokeExtensionMethod(reference.method, args)
        }
        else {
            throw new MissingMethodException("Channel.${name}", Object.class, args)
        }
    }

}
