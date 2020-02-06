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
    private Session session
    private boolean isLibrary
    private boolean anonymous

    @Deprecated
    IncludeDef( String module ) {
        this.modules = new ArrayList<>(1)
        this.modules << new Module(null, null)
        this.path = module
        this.anonymous = true
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

    IncludeDef setSession(Session session) {
        this.session = session
        return this
    }

    void define0(Map ownerParams) {
        checkValidPath(path)
        if( anonymous ) {
            // deprecated -- fallback on load
            log.warn "Anonymous module inclusion is deprecated -- Replace `include '${path}'` with `include { COMPONENT_NAME } from '${path}'`"
            loadModule(ownerParams)
        }
        else if( isLibrary ) {
            loadLibrary()
        }
        else {
            meta.defineModule(this)
        }
    }

    void loadLibrary() {
        // -- resolve the concrete against the current script
        final libraryFile = realModulePath(path)
        // -- load the library
        final libraryScript = loadLibrary0(libraryFile, session)
        // -- add it to the inclusions
        for( Module module : modules ) {
            meta.addModule(libraryScript, module.name, module.alias)
        }
    }

    /*
     * Note: this method invocation is injected during the Nextflow AST manipulation.
     * Do not use it explicitly.
     *
     * @param ownerParams The params in the owner context
     */
    void loadModule(Map ownerParams) {
        // -- resolve the concrete against the current script
        final moduleFile = realModulePath(path)
        // -- use the module specific params or default to the owner one if not provided
        final p = params != null ? params : ownerParams
        // -- load the module
        final moduleScript = loadModule0(moduleFile, p, session)
        // -- add it to the inclusions
        for( Module module : modules ) {
            meta.addModule(moduleScript, module.name, module.alias)
        }
    }

    @PackageScope
    ScriptMeta getMeta() { ScriptMeta.current() }

    @PackageScope
    Path getOwnerPath() { getMeta().getScriptPath() }

    @PackageScope
    boolean isLibrary() { isLibrary }

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
    @Memoized
    static BaseScript loadLibrary0(Path path, Session session) {
        // the execution of a library file has as side effect the registration of declared processes
        new ScriptParser(session)
                .setModule(true)
                .setBinding(new ScriptBinding())
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
            throw new IllegalModulePath("Missing include 'from' attribute")

        if( path instanceof Path && path.scheme != 'file' )
            throw new IllegalModulePath("Remote modules are not allowed -- Offending module: ${path.toUriString()}")

        final str = path.toString()
        if( !str.startsWith('/') && !str.startsWith('./') && !str.startsWith('../') )
            throw new IllegalModulePath("Module path must start with / or ./ prefix -- Offending module: $str")

        isLibrary = str.endsWith('.groovy') || str.endsWith('.nfl')
    }


}
