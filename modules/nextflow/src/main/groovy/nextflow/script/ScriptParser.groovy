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

import java.nio.file.Path

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.ast.NextflowXform
import nextflow.ast.OpXform
import nextflow.exception.ScriptCompilationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.io.ValueObject
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
/**
 * Parse a nextflow script class applied the required AST transformations
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ScriptParser {

    private ClassLoader classLoader

    private Session session

    private boolean module

    private Path scriptPath

    private BaseScript script

    private Object result

    private ScriptBinding binding

    private CompilerConfiguration config

    private String entryName

    ScriptParser(Session session) {
        this.session = session
        this.classLoader = session.classLoader
    }

    ScriptParser(ClassLoader loader) {
        this.classLoader = loader
    }

    ScriptParser setSession( Session session ) {
        this.session = session
        this.classLoader = session.classLoader
        return this
    }

    ScriptParser setModule(boolean value) {
        this.module = value
        return this
    }

    ScriptParser setEntryName(String name) {
        this.entryName = name
        return this
    }

    ScriptParser setBinding(ScriptBinding binding) {
        this.binding = binding
        return this
    }

    protected ClassLoader getClassLoader() { classLoader }

    protected Session getSession() { session }

    ScriptBinding getBinding() { binding }

    Object getResult() { result }

    BaseScript getScript() { script }

    CompilerConfiguration getConfig() {
        if( config )
            return config

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, groovy.transform.Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        importCustomizer.addImports( ValueObject.name )
        importCustomizer.addImport( 'channel', Channel.name )
        importCustomizer.addStaticStars( Nextflow.name )

        config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(OpXform))

        if( session?.debug )
            config.debug = true

        if( session?.classesDir )
            config.setTargetDirectory(session.classesDir.toFile())

        return config
    }

    /**
     * Creates a unique name for the main script class in order to avoid collision
     * with the implicit and user variables
     */
    protected String computeClassName(script) {
        final PREFIX = 'Script_'

        if( script instanceof Path ) {
            return FileHelper.getIdentifier(script,PREFIX)
        }

        if( script instanceof CharSequence ) {
            final hash = Hashing
                    .sipHash24()
                    .newHasher()
                    .putUnencodedChars(script.toString())
                    .hash()
            return PREFIX + hash.toString()
        }

        throw new IllegalArgumentException("Unknown script type: ${script?.getClass()?.getName()}")
    }

    private GroovyShell getInterpreter() {
        if( !binding && session )
            binding = session.binding
        if( !binding )
            throw new IllegalArgumentException("Missing Script binding object")

        return new GroovyShell(classLoader, binding, getConfig())
    }

    private ScriptParser parse0(String scriptText, Path scriptPath, GroovyShell interpreter) {
        this.scriptPath = scriptPath
        final String className = computeClassName(scriptText)
        try {
            final parsed = scriptPath && session.debug
                    ? interpreter.parse(scriptPath.toFile())
                    : interpreter.parse(scriptText, className)
            if( parsed !instanceof BaseScript ){
               throw new CompilationFailedException(0, null)
            }
            script = (BaseScript)parsed
            final meta = ScriptMeta.get(script)
            meta.setScriptPath(scriptPath)
            meta.setModule(module)
            meta.validate()
            return this
        }
        catch (CompilationFailedException e) {
            String type = module ? "Module" : "Script"
            String header = "$type compilation error\n- file : ${FilesEx.toUriString(scriptPath)}"
            String msg = e.message ?: header
            msg = msg != 'startup failed' ? msg : header
            msg = msg.replaceAll(/startup failed:\n/,'')
            msg = msg.replaceAll(~/$className(: \d+:\b*)?/, header+'\n- cause:')
            if( msg.contains "Unexpected input: '{'" ) {
                msg += "\nNOTE: If this is the beginning of a process or workflow, there may be a syntax error in the body, such as a missing or extra comma, for which a more specific error message could not be produced."
            }
            throw new ScriptCompilationException(msg, e)
        }
    }

    ScriptParser parse(String scriptText) {
        return parse0(scriptText, null, getInterpreter())
    }

    ScriptParser parse(Path scriptPath) {
        try {
            parse0(scriptPath.text, scriptPath, getInterpreter())
        }
        catch (IOException e) {
            throw new ScriptCompilationException("Unable to read script: '$scriptPath' -- cause: $e.message", e)
        }
        return this
    }

    ScriptParser runScript(String scriptText) {
        parse(scriptText)
        runScript()
        return this
    }

    ScriptParser runScript(Path scriptPath) {
        parse(scriptPath)
        runScript()
        return this
    }

    private void setupContext() {
        assert session
        binding.setSession(session)
        binding.setScriptPath(scriptPath)
        binding.setEntryName(entryName)
    }

    ScriptParser runScript() {
        assert script
        setupContext()
        result = script.run()
        return this
    }

}
