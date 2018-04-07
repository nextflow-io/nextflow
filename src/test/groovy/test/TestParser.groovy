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

package test
import java.nio.file.Path

import groovy.transform.InheritConstructors
import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.executor.Executor
import nextflow.processor.ProcessConfig
import nextflow.processor.ProcessFactory
import nextflow.processor.TaskProcessor
import nextflow.script.BaseScript
import nextflow.script.TaskBody
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

    Session session

    TestParser( Map config = null ) {
        session = config ? new Session(config) : new Session()
    }

    TestParser( Session session1 ) {
        session = session1
    }

    BaseScript parseScript ( String scriptText, Map map = null ) {
        parseScript(scriptText, new Binding(map))
    }

    BaseScript parseScript( String scriptText, Binding binding ) {

        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, groovy.transform.Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addStaticStars( Nextflow.name )

        def config = new CompilerConfiguration()
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers(importCustomizer)

        // extend the class-loader if required
        def gcl = new GroovyClassLoader()

        // run and wait for termination
        def groovy = new GroovyShell(gcl, binding, config)
        def script = groovy.parse( scriptText ) as BaseScript
        // initialize it
        script.setSession(session)
        script.setProcessFactory(new MockProcessFactory(script, session))
        // return it
        return script
    }

    TaskProcessor parseAndGetProcess( String scriptText ) {
        def script = parseScript(scriptText)
        script.run()
        return script.getTaskProcessor()
    }


    static TaskProcessor parseAndReturnProcess( String scriptText, Map map = [:] ) {
        def script = new TestParser().parseScript(scriptText, map)
        script.run()
        return script.getTaskProcessor()
    }


    @InheritConstructors
    static class MockProcessFactory extends ProcessFactory  {

        TaskProcessor newTaskProcessor( String name, Executor executor, Session session, BaseScript script, ProcessConfig config, TaskBody taskBody ) {
            new MockTaskProcessor(name, executor, session, script, config, taskBody)
        }

    }

    @InheritConstructors
    static class MockTaskProcessor extends TaskProcessor {
        @Override
        def run () { }
    }

}
