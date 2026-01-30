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
 */

package nextflow.script

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.NF
import nextflow.Session
import nextflow.config.ModulesConfig
import nextflow.config.RegistryConfig
import nextflow.exception.IllegalModulePath
import nextflow.exception.ScriptCompilationException
import nextflow.module.ModuleReference
import nextflow.module.ModuleResolver
import nextflow.pipeline.PipelineSpec
import nextflow.plugin.Plugins
import nextflow.plugin.extension.PluginExtensionProvider
import nextflow.script.parser.v1.ScriptLoaderV1
/**
 * Implements a script inclusion
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@EqualsAndHashCode
class IncludeDef {

    static final private String PLUGIN_PREFIX = 'plugin/'

    @Canonical
    static class Module {
        String name
        String alias
    }

    @PackageScope path
    @PackageScope List<Module> modules
    @PackageScope Map params
    @PackageScope Map addedParams
    private Session session

    protected IncludeDef(List<Module> modules) {
        this.modules = new ArrayList<>(modules)
    }

    IncludeDef from(Object path) {
        this.path = path
        return this
    }

    IncludeDef params(Map args) {
        log.warn1 "Include with `params()` is deprecated -- pass params as a workflow or process input instead"
        this.params = args != null ? new HashMap(args) : null
        return this
    }

    IncludeDef addParams(Map args) {
        log.warn1 "Include with `addParams()` is deprecated -- pass params as a workflow or process input instead"
        this.addedParams = args
        return this
    }

    IncludeDef setSession(Session session) {
        this.session = session
        return this
    }

    /**
     * Used internally by the script DSL to include modules
     * into a script.
     *
     * @param ownerParams The params in the including script context
     */
    void load0(ScriptBinding.ParamsMap ownerParams) {
        checkValidPath(path)
        if( path.toString().startsWith(PLUGIN_PREFIX) ) {
            loadPlugin0(path.toString().substring(PLUGIN_PREFIX.length()))
            return
        }
        // -- resolve the concrete against the current script
        final moduleFile = realModulePath(path).normalize()
        // -- load the module
        final moduleScript = NF.isSyntaxParserV2()
            ? loadModuleV2(moduleFile, ownerParams, session)
            : loadModuleV1(moduleFile, resolveParams(ownerParams), session)
        // -- add it to the inclusions
        for( Module module : modules ) {
            meta.addModule(moduleScript, module.name, module.alias)
        }
    }

    private Map resolveParams(ScriptBinding.ParamsMap ownerParams) {
        if( params!=null && addedParams!=null )
            throw new IllegalArgumentException("Include 'params' and 'addParams' option conflict -- check module: $path")
        if( params!=null )
            return params

        addedParams ? ownerParams.copyWith(addedParams) : ownerParams
    }

    @PackageScope
    ScriptMeta getMeta() { ScriptMeta.current() }

    @PackageScope
    Path getOwnerPath() { getMeta().getScriptPath() }

    /**
     * When using the strict syntax, the included script will already
     * have been compiled, so simply execute it to load its definitions.
     *
     * @param path    The included script path
     * @param params  The params of the including script
     * @param session The current workflow run
     */
    @PackageScope
    @Memoized
    static BaseScript loadModuleV2(Path path, Map params, Session session) {
        final script = ScriptMeta.getScriptByPath(path)
        if( !script )
            throw new IllegalStateException("Unable to find module script for path: $path")
        script.getBinding().setParams(params)
        script.run()
        return script
    }

    /**
     * When not using the strict syntax, compile and execute the
     * included script to load its definitions.
     *
     * @param path    The included script path
     * @param params  The params of the including script
     * @param session The current workflow run
     */
    @PackageScope
    @Memoized
    static BaseScript loadModuleV1(Path path, Map params, Session session) {
        final binding = new ScriptBinding() .setParams(params)
        new ScriptLoaderV1(session)
                .setModule(true)
                .setBinding(binding)
                .runScript(path)
                .getScript()
    }

    @PackageScope
    Path resolveModulePath(include) {
        assert include
        if( include.toString().startsWith('@') ) {
            return resolveRemoteModulePath(include.toString())
        }
        final result = include as Path
        if( result.isAbsolute() ) {
            if( result.scheme == 'file' ) return result
            throw new IllegalModulePath("Cannot resolve module path: ${result.toUriString()}")
        }

        return getOwnerPath().resolveSibling(include.toString())
    }

    @PackageScope
    Path resolveRemoteModulePath(String moduleName) {
        //TODO: Decide final location of modules currently in nextflow_spec.json.
        // Alternative: Use nextflow config. It requires to implement nextflow.config updater features
        // def modulesConfig = session.config.navigate('modules') as ModulesConfig
        final modulesConfig = getModuleConfig()
        final registryConfig = session.config.navigate('registry') as RegistryConfig
        // Create module resolver
        def resolver = new ModuleResolver(session.baseDir, modulesConfig, registryConfig)
        try {
            log.debug "Resolving remote module: ${moduleName}"

            // Parse module reference
            def reference = ModuleReference.parse(moduleName)

            // Resolve module (will auto-install if missing or version mismatch)
            def mainFile = resolver.resolve(reference, null,true)

            log.info "Module ${reference.nameWithoutPrefix} resolved to ${mainFile}"
            return mainFile
        } catch (Exception e) {
            throw new IllegalModulePath(
                "Failed to resolve remote module ${moduleName}: ${e.message}",
                e
            )
        }
    }

    @Memoized
    private ModulesConfig getModuleConfig(){
        def specFile = new PipelineSpec(session.getBaseDir())

        if( !specFile.exists() ) {
            log.warn1( "Remote module specified and 'nextflow_spec.json' not found." )
            return new ModulesConfig()
        }

        def modules = specFile.getModules()
        if (!modules || modules.isEmpty()) {
            log.warn1("Remote module specified and no modules configured in 'nextflow_spec.json'")
            return new ModulesConfig()
        }
        return new ModulesConfig(modules)
    }
    @PackageScope
    Path realModulePath(include) {
        def module = resolveModulePath(include)

        // check if exists a file with `.nf` extension
        if( !module.name.endsWith('.nf') ) {
            final extendedName = module.resolveSibling( "${module.name}.nf" )
            if( extendedName.exists() )
                return extendedName
        }
        if( module.isDirectory() ) {
            final target = module.resolve('main.nf')
            if( target.exists() ) {
                return target
            }
            throw new ScriptCompilationException("Include '$include' does not provide any module script -- the following path should contain a 'main.nf' script: '$module'" )
        }

        // check the file exists
        if( module.exists() )
            return module

        throw new NoSuchFileException("Can't find a matching module file for include: $include")
    }

    @PackageScope
    void checkValidPath(path) {
        if( !path )
            throw new IllegalModulePath("Missing module path attribute")

        if( path instanceof Path && path.scheme != 'file' )
            throw new IllegalModulePath("Remote modules are not allowed -- Offending module: ${path.toUriString()}")

        final str = path.toString()
        if( !str.startsWith('/') && !str.startsWith('./') && !str.startsWith('../') && !str.startsWith('plugin/') && !str.startsWith('@') )
            throw new IllegalModulePath("Module path must start with '/' , './' or '@' prefix -- Offending module: $str")

    }

    @PackageScope
    void loadPlugin0(String pluginId){
        if( pluginId.startsWith('/') )
            throw new IllegalArgumentException("Plugin Id in the 'include' declaration cannot start with a slash character - offending value: '$pluginId'")
        if( pluginId.contains('@') )
            throw new IllegalArgumentException("Plugin Id in the 'include' declaration cannot contain a specific version requirement - offending value: '$pluginId'")
        Plugins.startIfMissing(pluginId)
        if( !Plugins.isStarted(pluginId) )
            throw new IllegalArgumentException("Unable start plugin with Id '$pluginId'")
        final Map<String,String> declaredNames = this.modules.collectEntries {[it.name, it.alias ?: it.name]}
        log.debug "Loading included plugin extensions with names: $declaredNames; plugin Id: $pluginId"
        PluginExtensionProvider.INSTANCE().loadPluginExtensionMethods(pluginId, declaredNames)
    }

}
