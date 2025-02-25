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

import nextflow.exception.ScriptCompilationException
import nextflow.plugin.extension.PluginExtensionProvider
import nextflow.plugin.Plugins

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
import nextflow.exception.IllegalModulePath
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

    IncludeDef(TokenVar token, String alias=null) {
        def component = token.name; if(alias) component += " as $alias"
        def msg = "Unwrapped module inclusion is deprecated -- Replace `include $component from './MODULE/PATH'` with `include { $component } from './MODULE/PATH'`"
        if( NF.isDsl2() )
            throw new DeprecationException(msg)
        log.warn msg

        this.modules = new ArrayList<>(1)
        this.modules << new Module(token.name, alias)
    }

    protected IncludeDef(List<Module> modules) {
        this.modules = new ArrayList<>(modules)
    }

    /** only for testing purpose -- do not use */
    protected IncludeDef() { }

    IncludeDef from(Object path) {
        this.path = path
        return this
    }

    IncludeDef params(Map args) {
        log.warn "Include with `params()` is deprecated -- pass params as a workflow or process input instead"
        this.params = args != null ? new HashMap(args) : null
        return this
    }

    IncludeDef addParams(Map args) {
        log.warn "Include with `addParams()` is deprecated -- pass params as a workflow or process input instead"
        this.addedParams = args
        return this
    }

    IncludeDef setSession(Session session) {
        this.session = session
        return this
    }

    /*
     * Note: this method invocation is injected during the Nextflow AST manipulation.
     * Do not use it explicitly.
     *
     * @param ownerParams The params in the owner context
     */
    void load0(ScriptBinding.ParamsMap ownerParams) {
        checkValidPath(path)
        if( path.toString().startsWith(PLUGIN_PREFIX) ) {
            loadPlugin0(path.toString().substring(PLUGIN_PREFIX.length()))
            return
        }
        // -- resolve the concrete against the current script
        final moduleFile = realModulePath(path)
        // -- load the module
        final moduleScript = loadModule0(moduleFile, resolveParams(ownerParams), session)
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

    @PackageScope
    @Memoized
    static BaseScript loadModule0(Path path, Map params, Session session) {
        final binding = new ScriptBinding() .setParams(params)

        // the execution of a library file has as side effect the registration of declared processes
        new ScriptParser(session)
                .setModule(true)
                .setBinding(binding)
                .runScript(path)
                .getScript()
    }

    @PackageScope
    Path resolveModulePath(include) {
        assert include

        final result = include as Path
        if( result.isAbsolute() ) {
            if( result.scheme == 'file' ) return result
            throw new IllegalModulePath("Cannot resolve module path: ${result.toUriString()}")
        }

        return getOwnerPath().resolveSibling(include.toString())
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
        if( !str.startsWith('/') && !str.startsWith('./') && !str.startsWith('../') && !str.startsWith('plugin/') )
            throw new IllegalModulePath("Module path must start with / or ./ prefix -- Offending module: $str")

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
