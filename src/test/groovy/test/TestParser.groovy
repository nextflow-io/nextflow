/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package test
import java.nio.file.Path

import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session

import nextflow.ast.NextflowDSL
import nextflow.script.BaseScript
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
/**
 * An helper class to parse nextflow script snippets
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestParser {

    def session


    TestParser( Map session = null ) {
        this.session = session ? new Session(session) : new Session()
    }

    TestParser( Session session1 ) {
        this.session = session1
    }


    def parseScript ( String scriptText, Map map = null ) {
        parseScript(scriptText, new Binding(map))
    }

    def parseScript( String scriptText, Binding binding ) {

        binding.setVariable('__$TEST', true)
        binding.setVariable('__$session', session)

        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, groovy.transform.Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addStaticStars( Nextflow.name )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'splitter' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'operator' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'prioritySelector' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'select' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'selector' )

        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers(importCustomizer)
//        def preProcessor = new SourcePreProcessor()
//        preProcessor.lexerFactory = lexerFactory
//        config.setPluginFactory(preProcessor)

        // extend the class-loader if required
        def gcl = new GroovyClassLoader()

        // run and wait for termination
        def groovy = new GroovyShell(gcl, binding, config)
        groovy.parse( scriptText ) as BaseScript

    }


    static parse ( String scriptText, Map map = [:] ) {
        new TestParser().parseScript(scriptText, map)
    }

}
