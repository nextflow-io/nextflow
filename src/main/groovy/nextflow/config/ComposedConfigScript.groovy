/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
abstract class ComposedConfigScript extends Script {

    private static final Logger log = LoggerFactory.getLogger(ComposedConfigScript)

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
            includePath = owner.resolveSibling(includePath)
        }

        if( !includePath.exists() )
            throw new NoSuchFileException("Config file does not exist: ${includePath}")

        // -- set the required base script
        def config = new CompilerConfiguration()
        config.scriptBaseClass = ComposedConfigScript.class.name
        if( renderClosureAsString )
            config.addCompilationCustomizers(new ASTTransformationCustomizer(ConfigTransform))

        // -- setup the grengine instance
        def engine = new Grengine(this.class.classLoader,config)
        def clazz = engine.load(includePath.toUri().toURL())

        // -- push this file on the stack
        this.configStack.push(includePath)

        /*
         * here it is the magic code
         */
        def script = (ComposedConfigScript)clazz.newInstance()
        script.setConfigStack(this.configStack)
        script.setBinding(this.getBinding())
        script.metaClass.getProperty = { String name -> this.metaClass.getProperty(this, name) }
        script.metaClass.invokeMethod = { String name, args -> this.metaClass.invokeMethod(this, name, args)}
        script.&run.call()

        // remove the path from the stack
        this.configStack.pop()
    }

}
