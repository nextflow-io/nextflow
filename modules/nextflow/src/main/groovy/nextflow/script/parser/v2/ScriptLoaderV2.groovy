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

package nextflow.script.parser.v2

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.exception.ScriptCompilationException
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptLoader
import nextflow.script.ScriptMeta
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.SourceUnit
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Script parser/loader that uses the strict syntax.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ScriptLoaderV2 implements ScriptLoader {

    private Session session

    private BaseScript mainScript

    private boolean skipEntryFlow

    private Object result

    ScriptLoaderV2(Session session) {
        this.session = session
    }

    @Override
    ScriptLoaderV2 setEntryName(String name) {
        if( name )
            throw new IllegalArgumentException("The `-entry` option is not supported with the strict syntax -- use a param to run a named workflow from the entry workflow")
        return this
    }

    @Override
    ScriptLoaderV2 setModule(boolean value) {
        this.skipEntryFlow = value
        return this
    }

    @Override
    BaseScript getScript() {
        return mainScript
    }

    ScriptBinding getBinding() {
        mainScript.getBinding()
    }

    @Override
    Object getResult() {
        return result
    }

    @Override
    ScriptLoaderV2 parse(Path scriptPath) {
        try {
            parse0(scriptPath.text, scriptPath)
        }
        catch (IOException e) {
            throw new ScriptCompilationException("Unable to read script: '$scriptPath' -- cause: $e.message", e)
        }
        return this
    }

    ScriptLoaderV2 parse(String scriptText) {
        parse0(scriptText, null)
        return this
    }

    @Override
    ScriptLoaderV2 runScript() {
        assert session
        assert mainScript
        // capture the last statement of the snippet or entry workflow (used for testing)
        this.result = mainScript.run()
        return this
    }

    ScriptLoaderV2 runScript(String scriptText) {
        assert !mainScript
        parse(scriptText)
        runScript()
        return this
    }

    private void parse0(String scriptText, Path scriptPath) {
        final compiler = getCompiler()
        try {
            final compileResult = scriptPath
                ? compiler.compile(scriptPath.toFile())
                : compiler.compile(scriptText)

            this.mainScript = createScript(compileResult.main(), session.binding, scriptPath, skipEntryFlow)

            compileResult.modules().forEach((path, clazz) -> {
                createScript(clazz, new ScriptBinding(), path, true)
            })

            for( final name : compileResult.processNames() )
                ScriptMeta.addResolvedName(name)
        }
        catch( CompilationFailedException e ) {
            if( scriptPath )
                printErrors(scriptPath)
            throw new ScriptCompilationException("Script compilation failed", e)
        }
    }

    private void printErrors(Path path) {
        final errorListener = new StandardErrorListener('full', false)
        println()
        errorListener.beforeErrors()
        for( final message : compiler.getErrors() ) {
            final cause = message.getCause()
            final source = getSource(cause.getSourceLocator(), compiler)
            final filename = getRelativePath(source, path)
            errorListener.onError(cause, filename, source)
        }
        errorListener.afterErrors()
    }

    private SourceUnit getSource(String sourceLocator, ScriptCompiler compiler) {
        for( final su : compiler.getSources() ) {
            if( sourceLocator == su.getName() )
                return su
        }
        return null
    }

    private String getRelativePath(SourceUnit source, Path scriptPath) {
        final uri = source.getSource().getURI()
        return scriptPath.getParent().relativize(Path.of(uri)).toString()
    }

    private BaseScript createScript(Class clazz, ScriptBinding binding, Path path, boolean module) {
        final script = InvokerHelper.createScript(clazz, binding)
        if( script instanceof BaseScript ) {
            final meta = ScriptMeta.get(script)
            meta.setScriptPath(path)
            meta.setModule(module)
            binding.setScriptPath(path)
            binding.setSession(session)
            return script
        }
        throw new CompilationFailedException(0, null)
    }

    private ScriptCompiler compiler

    private ScriptCompiler getCompiler() {
        if( !compiler )
            compiler = new ScriptCompiler(session.debug, session.classesDir, session.getClassLoader())
        return compiler
    }

}
