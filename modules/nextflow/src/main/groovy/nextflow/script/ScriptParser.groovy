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

import java.nio.file.Path
import java.util.regex.Pattern

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.Channel
import nextflow.Nextflow
import nextflow.NextflowMeta
import nextflow.Session
import nextflow.ast.OpXform
import nextflow.ast.NextflowDSL
import nextflow.ast.NextflowXform
import nextflow.exception.ScriptCompilationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
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

    private static final Pattern DSL2_DECLARATION = ~/(?m)^\s*(nextflow\.(?:preview|enable)\.dsl\s*=\s*2)\s*;?\s*$/

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
        this.classLoader = session.getClassLoader()
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
        importCustomizer.addStaticStars( Nextflow.name )

        config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(OpXform))

        if( session && session.classesDir )
            config.setTargetDirectory(session.classesDir.toFile())

        return config
    }

    protected boolean isDsl2(String script) {
        script.find(DSL2_DECLARATION) != null
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
                    .murmur3_32()
                    .newHasher()
                    .putUnencodedChars(script.toString())
                    .hash()
            return PREFIX + hash.toString()
        }

        throw new IllegalArgumentException("Unknown script type: ${script?.getClass()?.getName()}")
    }

    GroovyShell getInterpreter() {
        if( !binding && session )
            binding = session.binding
        if( !binding )
            throw new IllegalArgumentException("Missing Script binding object")

        return new GroovyShell(classLoader, binding, getConfig())
    }

    ScriptParser parse(String scriptText, GroovyShell interpreter) {
        final String clazzName = computeClassName(scriptText)
        try {
            script = (BaseScript)interpreter.parse(scriptText, clazzName)
            final meta = ScriptMeta.get(script)
            meta.setScriptPath(scriptPath)
            meta.setModule(module)
            if( isDsl2(scriptText) )
                NextflowMeta.instance.enableDsl2()
            return this
        }
        catch (CompilationFailedException e) {
            String type = module ? "Module" : "Script"
            String header = "$type compilation error\n- file : ${FilesEx.toUriString(scriptPath)}"
            String msg = e.message ?: header
            msg = msg.replaceAll(/startup failed:\n/,'')
            msg = msg.replaceAll(~/$clazzName(: \d+:\b*)?/, header+'\n- cause:')
            throw new ScriptCompilationException(msg, e)
        }
    }


    ScriptParser parse(String scriptText) {
        def interpreter = getInterpreter()
        return parse(scriptText, interpreter)
    }

    ScriptParser parse(Path scriptPath) {
        this.scriptPath = scriptPath
        parse(scriptPath.text)
    }

    ScriptParser runScript(Path scriptPath) {
        this.scriptPath = scriptPath
        runScript(scriptPath.text)
        return this
    }

    ScriptParser runScript(String scriptText) {
        parse(scriptText)
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
