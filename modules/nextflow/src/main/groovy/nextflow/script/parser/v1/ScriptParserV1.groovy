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

package nextflow.script.parser.v1

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
import nextflow.script.BaseScript
import nextflow.script.ScriptParser
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
/**
 * Legacy script parser implement based on the Groovy parser.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ScriptParserV1 extends ScriptParser {

    private ClassLoader classLoader

    private CompilerConfiguration config

    ScriptParserV1(Session session) {
        super(session)
        this.classLoader = session.classLoader
    }

    @Override
    ScriptParserV1 setSession(Session session) {
        super.setSession(session)
        this.classLoader = session.classLoader
        return this
    }

    protected ClassLoader getClassLoader() { classLoader }

    @Override
    protected BaseScript parse0(String scriptText, Path scriptPath) {
        final interpreter = getInterpreter()
        final className = computeClassName(scriptText)
        try {
            final parsed = scriptPath && session.debug
                    ? interpreter.parse(scriptPath.toFile())
                    : interpreter.parse(scriptText, className)
            if( parsed !instanceof BaseScript )
               throw new CompilationFailedException(0, null)
            return (BaseScript)parsed
        }
        catch (CompilationFailedException e) {
            String type = isModule() ? "Module" : "Script"
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

    private GroovyShell getInterpreter() {
        if( !binding && session )
            binding = session.binding
        if( !binding )
            throw new IllegalArgumentException("Missing Script binding object")

        return new GroovyShell(classLoader, binding, getConfig())
    }

    CompilerConfiguration getConfig() {
        if( config )
            return config

        // define the imports
        final importCustomizer = new ImportCustomizer()
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
     * with implicit and user-defined variables.
     *
     * @param script
     */
    protected String computeClassName(script) {
        final String PREFIX = 'Script_'

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

}
