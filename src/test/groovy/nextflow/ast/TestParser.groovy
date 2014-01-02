/*
 * Copyright (c) 2012, the authors.
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

package nextflow.ast

import nextflow.Session
import nextflow.script.BaseScript
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestParser {

    def session = new Session()

    def Closure<CustomLexer> lexerFactory

    def parse ( String scriptText, Map map = [:] ) {
        parse(scriptText, new Binding(map))
    }

    def parse ( String scriptText, Binding binding ) {

        binding.setVariable('__$TEST', true)
        binding.setVariable('__$session', session)

        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(ProcessDefTransform))
        def preProcessor = new SourcePreProcessor()
        preProcessor.lexerFactory = lexerFactory
        config.setPluginFactory(preProcessor)

        // extend the class-loader if required
        def gcl = new GroovyClassLoader()

        // run and wait for termination
        def groovy = new GroovyShell(gcl, binding, config)
        groovy.parse( scriptText ) as BaseScript

    }


}
