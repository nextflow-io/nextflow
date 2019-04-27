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

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.transform.PackageScope
import nextflow.Session
import nextflow.exception.ProcessException
/**
 * Implements a script inclusion
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode
class IncludeDef {

    @PackageScope path
    @PackageScope String alias
    @PackageScope String name
    @PackageScope Map params = new LinkedHashMap(10)
    private Session session

    IncludeDef( String module ) {
        this.path = module
    }

    IncludeDef(TokenVar name, String alias=null) {
        this.name = name.name
        this.alias = alias
    }

    protected IncludeDef() {}

    IncludeDef from(Object path) {
        this.path = path
        return this
    }

    IncludeDef params(Map args) {
        this.params.putAll(args)
        return this
    }

    IncludeDef setSession(Session session) {
        this.session = session
        return this
    }

    void load() {
        if( !path )
            throw new IllegalArgumentException("Missing module path attribute")

        // -- resolve the concrete against the current script
        final moduleFile = realModulePath(path)
        // -- load the module
        def moduleScript = loadModule0(moduleFile, params, session)
        // -- add it to the inclusions
        meta.addModule(moduleScript, name, alias)
    }

    @PackageScope
    ScriptMeta getMeta() { ScriptMeta.current() }

    @PackageScope
    Path getOwnerPath() { getMeta().getScriptPath() }

    @PackageScope
    @Memoized
    static BaseScript loadModule0(Path path, Map params, Session session) {
        try {
            final binding = new ScriptBinding() .setParams(params)

            // the execution of a library file has as side effect the registration of declared processes
            new ScriptParser(session)
                    .setModule(true)
                    .setBinding(binding)
                    .runScript(path)
                    .getScript()
        }
        catch( ProcessException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Unable to load module file: $path -- cause: ${e.cause?.message ?: e.message}", e)
        }
    }

    @PackageScope
    Path resolveModulePath(include) {
        assert include

        final result = include as Path
        if( result.isAbsolute() ) {
            if( result.scheme == 'file' ) return result
            throw new IllegalArgumentException("Cannot resolve module path: ${result.toUriString()}")
        }

        return getOwnerPath().resolveSibling(include.toString())
    }

    @PackageScope
    Path realModulePath(include) {
        def module = resolveModulePath(include)

        // check if exists a file with `.nf` extension
        if( !module.name.contains('.') ) {
            def extendedName = module.resolveSibling( "${module.simpleName}.nf" )
            if( extendedName.exists() )
                return extendedName
        }

        // check the file exists
        if( module.exists() )
            return module

        throw new NoSuchFileException("Can't find a matching module file for include: $include")
    }

}
