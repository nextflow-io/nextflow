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

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.ast.NextflowXform
import nextflow.config.ConfigParser
import nextflow.config.StripSecretsXform
import nextflow.config.parser.ConfigParserPluginFactory
import nextflow.extension.Bolts
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer

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

    private Set<String> parsedProfiles

    private GroovyShell groovyShell

    @Override
    ConfigParserV2 setProfiles(List<String> profiles) {
        this.appliedProfiles = profiles
        return this
    }

    @Override
    Set<String> getProfiles() {
        return parsedProfiles
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

    /**
     * Parse the given script as a string and return the configuration object
     *
     * @param text
     * @param location
     */
    @Override
    ConfigObject parse(String text) {
        parse(text, null)
    }

    ConfigObject parse(String text, Path location) {
        final groovyShell = getGroovyShell()
        final script = (ConfigDsl) groovyShell.parse(text, uniqueClassName(text))
        if( location )
            script.setConfigPath(location)
        script.setIgnoreIncludes(ignoreIncludes)
        script.setRenderClosureAsString(renderClosureAsString)
        if( location )
            script.setConfigPath(location)

        script.setBinding(new Binding(bindingVars))
        script.setParams(paramVars)
        script.setProfiles(appliedProfiles)
        script.run()

        final target = script.getTarget()
        if( !target.params )
            target.remove('params')
        parsedProfiles = script.getParsedProfiles()
        return Bolts.toConfigObject(target)
    }

    @Override
    ConfigObject parse(File file) {
        return parse(file.toPath())
    }

    @Override
    ConfigObject parse(Path path) {
        return parse(path.text, path)
    }

    private GroovyShell getGroovyShell() {
        if( groovyShell )
            return groovyShell
        final classLoader = new GroovyClassLoader()
        final config = new CompilerConfiguration()
        config.setScriptBaseClass(ConfigDsl.class.getName())
        config.setPluginFactory(new ConfigParserPluginFactory())
        config.addCompilationCustomizers(new ASTTransformationCustomizer(ConfigToGroovyXform))
        if( stripSecrets )
            config.addCompilationCustomizers(new ASTTransformationCustomizer(StripSecretsXform))
        if( renderClosureAsString )
            config.addCompilationCustomizers(new ASTTransformationCustomizer(ClosureToStringXform))
        config.addCompilationCustomizers(new ASTTransformationCustomizer(NextflowXform))
        final importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        config.addCompilationCustomizers(importCustomizer)
        return groovyShell = new GroovyShell(classLoader, new Binding(), config)
    }

    /**
     * Creates a unique name for the config class in order to avoid collision
     * with config DSL
     *
     * @param text
     */
    private String uniqueClassName(String text) {
        def hash = Hashing
                .sipHash24()
                .newHasher()
                .putUnencodedChars(text)
                .hash()
        return "_nf_config_$hash"
    }

}
