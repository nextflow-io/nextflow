/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

    @Deprecated
    IncludeDef( String module ) {
        log.warn "Anonymous module inclusion is deprecated -- Replace `include '${module}'` with `include { MODULE_NAME } from '${module}'`"
        this.path = module
        this.modules = new ArrayList<>(1)
        this.modules << new Module(null,null)
    }

    IncludeDef(TokenVar name, String alias=null) {
        this.modules = new ArrayList<>(1)
        this.modules << new Module(name.name, alias)
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
        this.params = args != null ? new HashMap(args) : null
        return this
    }

    IncludeDef addParams(Map args) {
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
            def extendedName = module.resolveSibling( "${module.name}.nf" )
            if( extendedName.exists() )
                return extendedName
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
        if( !str.startsWith('/') && !str.startsWith('./') && !str.startsWith('../') )
            throw new IllegalModulePath("Module path must start with / or ./ prefix -- Offending module: $str")

    }


}
