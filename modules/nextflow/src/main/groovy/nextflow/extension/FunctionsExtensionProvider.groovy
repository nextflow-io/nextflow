package nextflow.extension

import nextflow.Global
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.script.FunctionDef

import java.lang.reflect.Modifier


/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
class FunctionsExtensionProvider {

    private static FunctionsExtensionProvider instance

    static FunctionsExtensionProvider INSTANCE() {
        if( instance != null )
            return instance
        return instance = new FunctionsExtensionProvider().install()
    }

    static void reset() {
        instance = null
    }

    FunctionsExtensionProvider install() {
        // configure as global instance
        return instance = this
    }

    /**
     * Load functions declared in the extension
     * @param pluginId
     * @param declaredNames
     * @return a list with methods not founded in the extension
     */
    List<FunctionDef> loadPluginFunctions(String pluginId, Map<String,String>declaredNames){
        final extensions= Plugins.getExtensionsInPluginId(FunctionExtensionPoint, pluginId)
        // If there are not functions factory them all declaredNames are channel extension
        if( !extensions ) {
            return declaredNames.keySet()
        }
        if (extensions.size() > 1)
            throw new AbortOperationException("Plugin '$pluginId' implements more than one extension point: ${extensions.collect(it -> it.class.getSimpleName()).join(',')}")
        FunctionExtensionPoint ext = extensions.first() as FunctionExtensionPoint
        def functionsLoaded = loadPluginFunctionMethods(ext, declaredNames)
        if( functionsLoaded.size() ){
            ext.checkInit(Global.session)
        }
        functionsLoaded
    }

    protected List<FunctionDef> loadPluginFunctionMethods( Object ext, Map<String,String> declaredNames){
        // find all function factories defined in the plugin
        final definedFunctions= getDeclaredFactoryExtensionMethods0(ext.getClass())
        List<FunctionDef> functionsLoaded = []
        declaredNames.each{ entry ->
            String realName = entry.key
            String aliasName = entry.value

            if( definedFunctions.contains(realName) ) {
                FunctionDef pluginFunctionDef = new FunctionDef(ext,realName, aliasName)
                functionsLoaded.add(pluginFunctionDef)
            }
        }
        functionsLoaded
    }

    static private Set<String> getDeclaredFactoryExtensionMethods0(Class clazz) {
        def result = new HashSet<String>(30)
        def methods = clazz.getDeclaredMethods()
        for( def handle : methods ) {
            // skip non-public methods
            if( !Modifier.isPublic(handle.getModifiers()) ) continue
            // skip static methods
            if( Modifier.isStatic(handle.getModifiers()) ) continue

            result.add(handle.name)
        }
        return result
    }
}
