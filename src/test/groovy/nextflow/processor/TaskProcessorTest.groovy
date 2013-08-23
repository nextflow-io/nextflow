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

package nextflow.processor
import nextflow.Session
import nextflow.executor.AbstractExecutor
import nextflow.executor.NopeExecutor
import nextflow.script.BaseScript
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskProcessorTest extends Specification {


    def testEnvironment() {

        setup:
        def home = File.createTempDir()
        def binFolder = new File(home,'bin')
        binFolder.mkdirs()

        when:
        def session = new Session([env: [X:"1", Y:"2"]])
        session.setBaseDir(home)
        def processor = new DummyProcessor(new NopeExecutor(), session, new DummyScript(), new TaskConfig(), {})
        def builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == binFolder.toString()

        when:
        session = new Session([env: [X:"1", Y:"2", PATH:'/some']])
        session.setBaseDir(home)
        processor = new DummyProcessor(new NopeExecutor(), session, new DummyScript(), new TaskConfig(), {})
        builder = new ProcessBuilder()
        builder.environment().putAll( processor.getProcessEnvironment() )

        then:
        noExceptionThrown()
        builder.environment().X == '1'
        builder.environment().Y == '2'
        builder.environment().PATH == "${binFolder.toString()}:/some"


        cleanup:
        home.deleteDir()

    }


    static class DummyProcessor extends TaskProcessor {

        DummyProcessor(AbstractExecutor executor, Session session, BaseScript script, TaskConfig taskConfig, Closure taskBlock) {
            super(executor, session, script, taskConfig, taskBlock)
        }

        @Override
        protected void createOperator() {
        }
    }

    static class DummyScript extends BaseScript {
        @Override
        Object run() {
            return null
        }
    }

}
