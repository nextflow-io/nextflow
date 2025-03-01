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

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.exception.ScriptCompilationException
import nextflow.extension.FilesEx
import nextflow.file.FileHelper
import nextflow.script.BaseScript
import nextflow.script.ScriptParser
import org.codehaus.groovy.control.CompilationFailedException
/**
 * Script parser implementation based on the Nextflow formal grammar.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ScriptParserV2 extends ScriptParser {

    private ScriptCompiler compiler

    ScriptParserV2(Session session) {
        super(session)
    }

    @Override
    protected BaseScript parse0(String scriptText, Path scriptPath) {
        final compiler = getCompiler()
        final className = computeClassName(scriptText)
        try {
            final parsed = scriptPath && session.debug
                    ? compiler.compile(scriptPath.toFile())
                    : compiler.compile(scriptText, className)
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
            throw new ScriptCompilationException(msg, e)
        }
    }

    private ScriptCompiler getCompiler() {
        if( !compiler ) {
            if( !binding && session )
                binding = session.binding
            if( !binding )
                throw new IllegalArgumentException("Missing Script binding object")
            compiler = new ScriptCompiler(session, binding)
        }

        return compiler
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
