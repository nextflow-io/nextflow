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

import groovy.transform.InheritConstructors
import nextflow.Session
import nextflow.executor.Executor
import nextflow.script.ProcessConfig
import nextflow.script.ProcessFactory
import nextflow.processor.TaskProcessor
import nextflow.script.BaseScript
import nextflow.script.ScriptBinding
import nextflow.script.ScriptParser
import nextflow.script.BodyDef
/**
 * An helper class to parse nextflow script snippets
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TestParser {

    Session session

    TestParser( Map config = null ) {
        session = config ? new TestSession(config) : new TestSession()
    }

    TestParser( Session session1 ) {
        session = session1
    }


    TaskProcessor parseAndGetProcess( String scriptText, Map binding=null ) {
        if( binding != null ) {
            session.binding = new ScriptBinding(binding)
        }

        session.init(null,null)
        
        return new ScriptParser(session)
                .runScript(scriptText)
                .getScript()
                .getTaskProcessor()
    }


    static TaskProcessor parseAndReturnProcess( String scriptText, Map map = null ) {
        new TestParser().parseAndGetProcess(scriptText, map)
    }


    static class MockProcessFactory extends ProcessFactory {

        BaseScript script
        Session session

        MockProcessFactory(BaseScript script, Session session) {
            super(script,session)
            this.script = script
            this.session = session
        }

        @Override
        TaskProcessor newTaskProcessor(String name, Executor executor, ProcessConfig config, BodyDef taskBody ) {
            new MockTaskProcessor(name, executor, session, script, config, taskBody)
        }

    }

    @InheritConstructors
    static class MockTaskProcessor extends TaskProcessor {
        @Override
        def run () { }
    }

    @InheritConstructors
    static class TestSession extends Session {

        @Override
        ProcessFactory newProcessFactory(BaseScript script) {
            return new MockProcessFactory(script, this)
        }
    }
}
