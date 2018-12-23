/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import java.nio.file.NoSuchFileException
import java.nio.file.Path

import ch.grengine.Grengine
import nextflow.file.FileHelper
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
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
     * Implements the config file include
     */
    def includeConfig( includeFile ) {
        assert includeFile

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

        def configText
        try {
            configText = includePath.getText()
        }
        catch( IOException e ) {
            if( !includePath.exists() )
                throw new NoSuchFileException("Config file does not exist: ${includePath.toUriString()}")
            else
                throw new IOException("Cannot read config file include: ${includePath.toUriString()}", e)
        }

        // -- set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ConfigBase.class.name
        def params = [:]
        if( renderClosureAsString )
            params.put('renderClosureAsString', true)
        config.addCompilationCustomizers(new ASTTransformationCustomizer(params, ConfigTransform))

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

}
