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

package nextflow.config.parser.v2

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.config.ConfigParser
import nextflow.exception.ConfigParseException
import nextflow.extension.Bolts
import nextflow.script.parser.v2.StandardErrorListener
import org.codehaus.groovy.control.CompilationFailedException
import org.codehaus.groovy.control.SourceUnit

/**
 * The parser for Nextflow config files.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ConfigParserV2 implements ConfigParser {

    private Map bindingVars = [:]

    private Map paramVars = [:]

    private boolean ignoreIncludes = false

    private boolean renderClosureAsString = false

    private boolean strict = true

    private boolean stripSecrets

    private List<String> appliedProfiles

    private Set<String> declaredProfiles

    private Map<String,Object> declaredParams

    private GroovyShell groovyShell

    @Override
    ConfigParserV2 setProfiles(List<String> profiles) {
        this.appliedProfiles = profiles
        return this
    }

    @Override
    ConfigParserV2 setIgnoreIncludes(boolean value) {
        this.ignoreIncludes = value
        return this
    }

    @Override
    ConfigParserV2 setRenderClosureAsString(boolean value) {
        this.renderClosureAsString = value
        return this
    }

    @Override
    ConfigParserV2 setStrict(boolean value) {
        this.strict = value
        return this
    }

    @Override
    ConfigParser setStripSecrets(boolean value) {
        this.stripSecrets = value
        return this
    }

    @Override
    ConfigParserV2 setBinding(Map vars) {
        this.bindingVars = vars
        return this
    }

    @Override
    ConfigParserV2 setParams(Map vars) {
        // deep clone the map to prevent side-effect
        // see https://github.com/nextflow-io/nextflow/issues/1923
        this.paramVars = Bolts.deepClone(vars)
        return this
    }

    @Override
    Set<String> getDeclaredProfiles() {
        return declaredProfiles
    }

    @Override
    Map<String,Object> getDeclaredParams() {
        return declaredParams
    }

    /**
     * Parse the given script as a string and return the configuration object
     *
     * @param text
     * @param path
     */
    @Override
    ConfigObject parse(String text) {
        parse(text, null)
    }

    ConfigObject parse(String text, Path path) {
        final compiler = getCompiler()
        try {
            final script = (ConfigDsl) compiler.compile(text, path)
            script.setBinding(new Binding(bindingVars))
            if( path )
                script.setConfigPath(path)
            script.setIgnoreIncludes(ignoreIncludes)
            script.setParams(paramVars)
            script.setProfiles(appliedProfiles)
            script.setRenderClosureAsString(renderClosureAsString)
            script.run()

            final target = script.getTarget()
            declaredProfiles = script.getDeclaredProfiles()
            declaredParams = script.getDeclaredParams()
            return Bolts.toConfigObject(target)
        }
        catch( CompilationFailedException e ) {
            if( path )
                printErrors(path)
            throw new ConfigParseException("Config parsing failed", e)
        }
    }

    private void printErrors(Path path) {
        final source = compiler.getSource()
        final errorListener = new StandardErrorListener('full', false)
        println()
        errorListener.beforeErrors()
        for( final message : compiler.getErrors() ) {
            final cause = message.getCause()
            final filename = getRelativePath(source, path)
            errorListener.onError(cause, filename, source)
        }
        errorListener.afterErrors()
    }

    private String getRelativePath(SourceUnit source, Path path) {
        final uri = source.getSource().getURI()
        return path.getParent().relativize(Path.of(uri)).toString()
    }

    @Override
    ConfigObject parse(File file) {
        return parse(file.toPath())
    }

    @Override
    ConfigObject parse(Path path) {
        return parse(path.text, path)
    }

    private ConfigCompiler compiler

    private ConfigCompiler getCompiler() {
        if( !compiler )
            compiler = new ConfigCompiler(renderClosureAsString, stripSecrets)
        return compiler
    }

}
