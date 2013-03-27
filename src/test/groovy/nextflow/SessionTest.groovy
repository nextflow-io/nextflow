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

package nextflow

import nextflow.processor.LocalScriptProcessor
import nextflow.processor.NopeScriptProcessor
import nextflow.processor.OgeScriptProcessor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SessionTest extends Specification {


    def 'test createProcessor - bindOnTermination attribute' () {

        setup:
        def session = new Session()
        session.processorClass = NopeScriptProcessor

        when:
        def processor1 = session.createProcessor()
        def processor2 = session.createProcessor(false)
        def processor3 = session.createProcessor(true)

        then:
        !processor1.getBindOnTermination()
        !processor2.getBindOnTermination()
        processor3.getBindOnTermination()
    }

    def 'test createProcessor with attributes' () {

        setup:
        def taskConf = [
                echo: true,
                shareWorkDir: true,
                shell: 'zsh',
                threads: 3,
                environment: [a:1, b:2,c:3],
                validExitCodes: [1,2,3]
        ]
        def session = new Session( [task: taskConf] )

        when:
        def p = session.createProcessor()

        then:
        p.echo
//        p.shareWorkDirectory
        p.threads == 3
        p.environment == [a:1, b:2,c:3]
        p.validExitCodes == [1,2,3]
        p.shell == 'zsh'
    }


    def 'test loadProcessorClass' () {

        expect:
        new Session().loadProcessorClass(null) == LocalScriptProcessor
        new Session().loadProcessorClass('local') == LocalScriptProcessor
        new Session().loadProcessorClass('sge') == OgeScriptProcessor
        new Session().loadProcessorClass( OgeScriptProcessor.name ) == OgeScriptProcessor

    }

    def 'test new Session with config' () {
        expect:
        new Session().processorClass == LocalScriptProcessor
        new Session([task: [processor:'oge']]).processorClass == OgeScriptProcessor
    }

//
//    def 'test tasksReport' () {
//
//        given:
//        def session = new Session()
//        def processor = new LocalScriptProcessor(session).name('Hola')
//
//        TaskDef task1 = new TaskDef(id: 1, status: TaskDef.Status.RUNNING)
//        TaskDef task2 = new TaskDef(id: 2, status: TaskDef.Status.TERMINATED, exitCode: 1)
//
//        session.tasks.put(processor, task1)
//        session.tasks.put(processor, task2)
//
//        when:
//        print session.tasksReport()
//
//        then:
//        noExceptionThrown()
//
//    }



}
