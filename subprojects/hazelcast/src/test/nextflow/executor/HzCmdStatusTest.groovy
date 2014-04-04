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

package nextflow.executor

import java.nio.file.Paths

import nextflow.processor.DelegateMap
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzCmdStatusTest extends Specification {

    def testFactoryMethods() {

        setup:
        def ctx = [x:1, y:2]
        def sessionId = UUID.randomUUID()
        def map = [ dehydrate:{}, getHolder: {ctx} ] as DelegateMap

        def closure = { -> }
        closure.delegate = map
        def task = new TaskRun(id:'123', workDirectory: Paths.get('xxx'), code: closure)
        def cmd = new HzCmdCall( sessionId, task )

        /*
         * *start* message
         */
        when:
        def message = HzCmdStatus.start( cmd, 'xxx-yyy' )
        then:
        message.sessionId == sessionId
        message.taskId == '123'
        message.memberId == 'xxx-yyy'
        message.event == HzCmdStatus.Event.START

        /*
         * *result* message
         */
        when:
        message = HzCmdStatus.result( cmd, 33, new RuntimeException() )
        then:
        message.sessionId == sessionId
        message.taskId == '123'
        message.value == 33
        message.error instanceof RuntimeException
        message.context == ctx
        message.event == HzCmdStatus.Event.COMPLETE

        /*
         * error message
         */
        when:
        message = HzCmdStatus.error( sessionId, 'task_x' )
        then:
        message.sessionId == sessionId
        message.taskId == 'task_x'
        message.event == HzCmdStatus.Event.COMPLETE

    }

    def testConstructor() {

        when:
        def uuid = UUID.randomUUID()
        def err = new RuntimeException()
        def cmd = new HzCmdStatus( sessionId: uuid, taskId: '2', value: '3', error: err, context: [a:1,b:2], memberId: 'abc' )
        then:
        cmd.sessionId == uuid
        cmd.taskId == '2'
        cmd.value == '3'
        cmd.error == err
        cmd.context == [a:1,b:2]
        cmd.memberId == 'abc'

        when:
        def id2 = UUID.randomUUID()
        def copy1 = cmd.copyWith( sessionId: id2 )
        then:
        copy1.sessionId == id2
        copy1.taskId == '2'
        copy1.value == '3'
        copy1.error == err
        copy1.context == [a:1,b:2]
        copy1.memberId == 'abc'
        copy1 != cmd

        when:
        def copy2 = cmd.copyWith( taskId: 'omega' )
        then:
        copy2.sessionId == uuid
        copy2.taskId == 'omega'
        copy2.value == '3'
        copy2.error == err
        copy2.context == [a:1,b:2]
        copy2.memberId == 'abc'
        copy2 != cmd

    }

}
