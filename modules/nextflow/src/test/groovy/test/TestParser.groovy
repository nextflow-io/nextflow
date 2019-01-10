/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package test
import java.nio.file.Path

import groovy.transform.InheritConstructors
import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.ast.NextflowXform
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
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))
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
