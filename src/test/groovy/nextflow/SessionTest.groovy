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

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SessionTest extends Specification {

    def 'test parseEnvFile' () {

        setup:
        def map = [X:'hola',Y:'hello']

        def file = File.createTempFile('test','env')
        file.deleteOnExit()
        file.text = '''
            A=1
            B=${X}_123
            C=${Y}
            D=${B}:${C}
            Q='$xyz'
            P="${Y}"
            Z=${missing_variable}
            '''
        when:
        def session = new Session()
        def result = session.parseEnvironmentFile(file, map)


        then:
        result['A'] == '1'
        result['B'] == 'hola_123'
        result['C'] == 'hello'
        result['D'] == "hola_123:hello"
        result['Q'] == '$xyz'
        result['P'] == 'hello'
        result['Z'] == ''



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
