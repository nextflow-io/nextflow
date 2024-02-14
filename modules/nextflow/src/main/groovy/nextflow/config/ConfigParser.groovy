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

package nextflow.config

import java.nio.file.Path

import ch.artecat.grengine.Grengine
import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.ast.NextflowXform
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
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ConfigParser {

    private Map bindingVars = [:]

    private Map paramVars = [:]

    private boolean ignoreIncludes = false

    private boolean renderClosureAsString = false

    private boolean strict = true

    private Grengine grengine

    /**
     * Returns the profile names defined in the config file
     *
     * @return The set of profile names.
     */
    Set<String> getProfileNames() { profileNames }

    private Grengine getGrengine() {
        if( grengine ) {
            return grengine
        }

        // set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ConfigDsl.class.name
        config.addCompilationCustomizers(new ASTTransformationCustomizer(ConfigTransform))
        config.addCompilationCustomizers(new ASTTransformationCustomizer(NextflowXform))
        //  add implicit types
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( ClosureWithSource.name )
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        config.addCompilationCustomizers(importCustomizer)
        return grengine = new Grengine(config)
    }

    /**
     * Disable parsing of {@code includeConfig} directive
     *
     * @param value A boolean value, when {@code true} includes are disabled
     * @return The {@link ConfigParser} object itself
     */
    ConfigParser setIgnoreIncludes(boolean value) {
        this.ignoreIncludes = value
        return this
    }

    ConfigParser setRenderClosureAsString(boolean value) {
        this.renderClosureAsString = value
        return this
    }

    ConfigParser setStrict(boolean value) {
        this.strict = value
        return this
    }

    ConfigParser setBinding(Map vars) {
        this.bindingVars = vars
        return this
    }

    ConfigParser setParams(Map vars) {
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
    ConfigObject parse(String text, Path location=null) {
        final grengine = getGrengine()
        final encoded = "'${text.bytes.encodeBase64().toString()}'"
        final dsl = (ConfigDsl)grengine.load(encoded, createUniqueName(text)).newInstance()
        dsl.setIgnoreIncludes(ignoreIncludes)
        dsl.setRenderClosureAsString(renderClosureAsString)
        dsl.setStrict(strict)
        if( location )
            dsl.setConfigPath(location)

        dsl.setBinding(new Binding(bindingVars))
        dsl.setParams(paramVars)
        dsl.run()

        return Bolts.toConfigObject(dsl.getTarget())
    }

    ConfigObject parse(File file) {
        return parse(file.toPath())
    }

    ConfigObject parse(Path path) {
        return parse(path.text, path)
    }

}
