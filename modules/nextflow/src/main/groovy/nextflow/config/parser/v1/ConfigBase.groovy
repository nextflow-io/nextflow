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

package nextflow.config.parser.v1

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import ch.artecat.grengine.Grengine
import groovy.transform.Memoized
import nextflow.SysEnv
import nextflow.config.StripSecretsXform
import nextflow.exception.IllegalConfigException
import nextflow.file.FileHelper
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Extends a {@link Script} class adding a method that allows a configuration
 * file to include other configuration files.
 * <p>
 * This class is used to base class when parsing the groovy config object
 * <p>
 *
 * Based on
 * http://naleid.com/blog/2009/07/30/modularizing-groovy-config-files-with-a-dash-of-meta-programming
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class ConfigBase extends Script {

    private static final Logger log = LoggerFactory.getLogger(ConfigBase)

    private Stack<Path> configStack

    private boolean ignoreIncludes

    private boolean renderClosureAsString

    private boolean stripSecrets

    protected void setStripSecrets( boolean value ) {
        this.stripSecrets = value
    }

    protected void setIgnoreIncludes( boolean value ) {
        this.ignoreIncludes = value
    }

    protected void setRenderClosureAsString( boolean value ) {
        this.renderClosureAsString = value
    }

    protected void setConfigPath(Path path) {
        if( configStack == null )
            configStack = new Stack<>()
        configStack.push(path)
    }

    protected void setConfigStack(Stack<Path> stack) {
        this.configStack = stack
    }

    /**
     * Get the value of an environment variable from the launch environment.
     *
     * @param name
     *      The environment variable name to be referenced
     * @return
     *      The value associate with the specified variable name or {@code null} if the variable does not exist.
     */
    String env(String name) {
        return SysEnv.get(name)
    }

    /**
     * Implements the config file include
     */
    def includeConfig( includeFile ) {
        if( !includeFile )
            throw new IllegalConfigException("includeConfig argument cannot be empty")

        if( ignoreIncludes )
            return

        if( configStack == null )
            configStack = new Stack<>()

        Path owner = configStack ? this.configStack.peek() : null
        Path includePath = FileHelper.asPath(includeFile.toString())
        log.trace "Include config file: $includeFile [parent: $owner]"

        if( !includePath.isAbsolute() && owner ) {
            includePath = owner.resolveSibling(includeFile.toString())
        }

        def configText = readConfigFile(includePath)

        // -- set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ConfigBase.class.name
        if( stripSecrets )
            config.addCompilationCustomizers(new ASTTransformationCustomizer(StripSecretsXform))
        def params = [:]
        if( renderClosureAsString )
            params.put('renderClosureAsString', true)
        config.addCompilationCustomizers(new ASTTransformationCustomizer(params, ConfigTransform))
        //  -- add implicit types
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( Duration.name )
        importCustomizer.addImports( MemoryUnit.name )
        config.addCompilationCustomizers(importCustomizer)

        // -- setup the grengine instance
        def engine = new Grengine(this.class.classLoader,config)
        def clazz = engine.load(configText)

        // -- push this file on the stack
        this.configStack.push(includePath)

        /*
         * here it is the magic code
         */
        def script = (ConfigBase)clazz.newInstance()
        script.setConfigStack(this.configStack)
        script.setBinding(this.getBinding())
        script.metaClass.getProperty = { String name -> this.metaClass.getProperty(this, name) }
        script.metaClass.invokeMethod = { String name, args -> this.metaClass.invokeMethod(this, name, args)}
        script.&run.call()

        // remove the path from the stack
        this.configStack.pop()
    }

    /**
     * Reads the content of a config file. The result is cached to
     * avoid multiple reads.
     *
     * @param includePath The config file path
     * @return The file content
     */
    @Memoized
    protected static String readConfigFile(Path includePath) {
        def configText
        try {
            configText = includePath.getText()
        }
        catch (NoSuchFileException | FileNotFoundException ignored) {
            throw new NoSuchFileException("Config file does not exist: ${includePath.toUriString()}")
        }
        catch (IOException e) {
            throw new IOException("Cannot read config file include: ${includePath.toUriString()}", e)
        }
        return configText
    }

}
