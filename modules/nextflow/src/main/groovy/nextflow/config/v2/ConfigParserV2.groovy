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

package nextflow.config.v2

import java.nio.file.Path

import ch.artecat.grengine.Grengine
import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.ast.NextflowXform
import nextflow.config.ConfigParser
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

    private List<String> appliedProfiles

    private Set<String> parsedProfiles

    private Grengine grengine

    @Override
    ConfigParserV2 setProfiles(List<String> profiles) {
        this.appliedProfiles = profiles
        return this
    }

    @Override
    Set<String> getProfiles() {
        return parsedProfiles
    }

    private Grengine getGrengine() {
        if( grengine ) {
            return grengine
        }

        // set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ConfigDsl.class.name
        config.setPluginFactory(new ConfigParserPluginFactory())
        config.addCompilationCustomizers(new ASTTransformationCustomizer(NextflowXform))
        //  add implicit types
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( ClosureWithSource.name )
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        config.addCompilationCustomizers(importCustomizer)
        return grengine = new Grengine(config)
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
     * Creates a unique name for the config class in order to avoid collision
     * with top level configuration scopes
     *
     * @param text
     */
    private String createUniqueName(String text) {
        def hash = Hashing
                .murmur3_32()
                .newHasher()
                .putUnencodedChars(text)
                .hash()
        return "_nf_config_$hash"
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
        final grengine = getGrengine()
        final dsl = (ConfigDsl)grengine.load(text, createUniqueName(text)).newInstance()
        dsl.setIgnoreIncludes(ignoreIncludes)
        dsl.setRenderClosureAsString(renderClosureAsString)
        dsl.setStrict(strict)
        if( location )
            dsl.setConfigPath(location)

        dsl.setBinding(new Binding(bindingVars))
        dsl.setParams(paramVars)
        dsl.run()

        final result = Bolts.toConfigObject(dsl.getTarget())
        final profiles = (result.profiles ?: [:]) as ConfigObject
        parsedProfiles = profiles.keySet()
        if( appliedProfiles ) {
            for( def profile : appliedProfiles )
                if( profile in parsedProfiles )
                    result.merge(profiles[profile] as ConfigObject)
            result.remove('profiles')
        }

        return result
    }

    @Override
    ConfigObject parse(File file) {
        return parse(file.toPath())
    }

    @Override
    ConfigObject parse(Path path) {
        return parse(path.text, path)
    }

}
