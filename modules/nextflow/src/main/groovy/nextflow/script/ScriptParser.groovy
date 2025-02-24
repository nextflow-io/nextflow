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

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.exception.ScriptCompilationException
/**
 * Interface for parsing and executing a Nextflow script.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
abstract class ScriptParser {

    private ClassLoader classLoader

    private Session session

    private boolean module

    private String entryName

    private ScriptBinding binding

    private Path scriptPath

    private BaseScript script

    private Object result

    ScriptParser(Session session) {
        this.session = session
        this.classLoader = session.classLoader
    }

    ScriptParser setSession(Session session) {
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

    protected boolean isModule() { module }

    ScriptBinding getBinding() { binding }

    Object getResult() { result }

    BaseScript getScript() { script }

    ScriptParser parse(String scriptText) {
        return parse(scriptText, null)
    }

    ScriptParser parse(Path scriptPath) {
        return parse(scriptPath.text, scriptPath)
    }

    protected ScriptParser parse(String scriptText, Path scriptPath) {
        try {
            this.scriptPath = scriptPath
            this.script = parse0(scriptText, scriptPath)
            final meta = ScriptMeta.get(script)
            meta.setScriptPath(scriptPath)
            meta.setModule(module)
            meta.validate()
        }
        catch (IOException e) {
            throw new ScriptCompilationException("Unable to read script: '$scriptPath' -- cause: $e.message", e)
        }
        return this
    }

    abstract protected BaseScript parse0(String scriptText, Path scriptPath)

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
