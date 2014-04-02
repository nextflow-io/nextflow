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
class HzCmdNotifyTest extends Specification {

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
        def message = HzCmdNotify.start( cmd, 'xxx-yyy' )
        then:
        message.sessionId == sessionId
        message.taskId == '123'
        message.memberId == 'xxx-yyy'
        message.event == HzCmdNotify.Event.START

        /*
         * *result* message
         */
        when:
        message = HzCmdNotify.result( cmd, 33, new RuntimeException() )
        then:
        message.sessionId == sessionId
        message.taskId == '123'
        message.value == 33
        message.error instanceof RuntimeException
        message.context == ctx
        message.event == HzCmdNotify.Event.COMPLETE

        /*
         * error message
         */
        when:
        message = HzCmdNotify.error( sessionId, 'task_x' )
        then:
        message.sessionId == sessionId
        message.taskId == 'task_x'
        message.event == HzCmdNotify.Event.COMPLETE


    }

}
