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

package nextflow.script

import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session
import nextflow.processor.NopeScriptProcessor
import org.codehaus.groovy.control.MultipleCompilationErrorsException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowRunnerTest extends Specification {

    def 'test task' () {
        setup:
        def session = new Session()
        session.processorClass = NopeScriptProcessor
        def runner = new CliRunner(session)

        /*
         * Test a task with a very simple body.
         * For testing purposes the processor just return the script itself as result
         */
        when:
        def script =
            """
            task('task1')  {
                "echo Hello world"
            }
            """

        def result = runner.execute(script)

        then:
        result == "echo Hello world"

    }


    def 'test task with assignment' () {
        setup:
        def session = new Session()
        session.processorClass = NopeScriptProcessor
        def runner = new CliRunner(session)

        when:
        def script =
            '''
            task('task2')  {

                input x: 1
                input y: [3,4]

                def z = 'Hello world!'

                """echo $z"""
            }

            '''
        def result = runner.execute(script)

        then:
        result == 'echo Hello world!'
        runner.getScript().taskProcessor().getName() == 'task2'
        runner.getScript().taskProcessor().getInput('x').getVal() == 1
        runner.getScript().taskProcessor().getInput('y') instanceof DataflowQueue
        runner.getScript().taskProcessor().getOutput('-') instanceof DataflowWriteChannel

    }

    def 'test task with syntax error' () {
        setup:
        def session = new Session()
        session.processorClass = NopeScriptProcessor
        def runner = new CliRunner(session)

        /*
         * this declaration returns a syntax error because the task code block
         * does not terminate with a script to execute
         */
        when:
        def script =
            """
            task('task1')  {
                input x: 'hola'
            }
            """
        runner.execute(script)

        then:
        thrown(MultipleCompilationErrorsException)
    }
}
