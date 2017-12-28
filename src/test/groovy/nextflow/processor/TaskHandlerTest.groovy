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

package nextflow.processor

import nextflow.Session
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskHandlerTest extends Specification {

    def static final long KB = 1024

    def 'test get trace record'() {

        given:
        def traceText =  '''
        pid state %cpu %mem vmem rss peak_vmem peak_rss rchar wchar syscr syscw read_bytes write_bytes
        1 0 10 20 11084 1220 21084 2220 4790 12 11 1 20 30
        '''
                .leftTrim()

        def folder = TestHelper.createInMemTempDir()
        folder.resolve( TaskRun.CMD_TRACE ).text = traceText

        def config = [
                tag: 'seq_x',
                container: 'ubuntu',
                queue: 'longjobs',
                cpus: 2,
                time: '1 hour',
                disk: '100 GB',
                memory: '4 GB'
        ]
        def task = new TaskRun(id: new TaskId(100), workDir: folder, name:'task1', exitStatus: 127, config: config  )
        task.metaClass.getHashLog = { "5d5d7ds" }
        task.processor = Mock(TaskProcessor)
        task.processor.getSession() >> new Session()
        task.processor.getName() >> 'TheProcessName'
        task.processor.getProcessEnvironment() >> [FOO:'hola', BAR: 'mundo']
        task.context = new TaskContext(Mock(Script), [:], 'none')

        def handler = [:] as TaskHandler
        handler.task = task
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500

        when:
        def trace = handler.getTraceRecord()

        then:
        trace.task_id == 100
        trace.status == 'COMPLETED'
        trace.hash == '5d5d7ds'
        trace.name == 'task1'
        trace.exit == 127
        trace.submit == 1000
        trace.start == 1500
        trace.process == 'TheProcessName'
        trace.tag == 'seq_x'
        trace.'%cpu' == 1.0f
        trace.'%mem' == 2.0f
        trace.rss == 1220 * KB
        trace.vmem == 11084 * KB
        trace.peak_rss == 2220 * KB
        trace.peak_vmem == 21084 * KB
        trace.rchar == 4790
        trace.wchar == 12
        trace.syscr == 11
        trace.syscw ==  1
        trace.read_bytes == 20
        trace.write_bytes == 30
        trace.queue == 'longjobs'
        trace.cpus == 2
        trace.time == Duration.of('1 hour').toMillis()
        trace.memory == MemoryUnit.of('4 GB').toBytes()
        trace.disk == MemoryUnit.of('100 GB').toBytes()
        trace.env == 'FOO=hola\nBAR=mundo\n'

        // check get method
        trace.getFmtStr('%cpu') == '1.0%'
        trace.getFmtStr('%mem') == '2.0%'
        trace.getFmtStr('rss') == '1.2 MB'
        trace.getFmtStr('vmem') == '10.8 MB'
        trace.getFmtStr('peak_rss') == '2.2 MB'
        trace.getFmtStr('peak_vmem') == '20.6 MB'
        trace.getFmtStr('time') == '1h'
        trace.getFmtStr('memory') == '4 GB'
        trace.getFmtStr('disk') == '100 GB'

        when:
        handler = [:] as TaskHandler
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500
        handler.task = task.clone()
        handler.task.failed = true
        then:
        handler.getTraceRecord().status == 'FAILED'


        when:
        handler = [:] as TaskHandler
        handler.status = TaskStatus.COMPLETED
        handler.submitTimeMillis = 1000
        handler.startTimeMillis = 1500
        handler.task = task.clone()
        handler.task.aborted = true
        then:
        handler.getTraceRecord().status == 'ABORTED'

    }


}
