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

package nextflow.processor

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

        when:
        def handler = [:] as TaskHandler
        handler.task = new TaskRun(workDir: folder)
        handler.status = TaskHandler.Status.COMPLETED
        def trace = handler.getTraceRecord()

        then:
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

        // check get method
        trace.get('%cpu', null) == '1.0%'
        trace.get('%mem', null) == '2.0%'
        trace.get('rss', null ) == '1.2 MB'
        trace.get('vmem', null ) == '10.8 MB'
        trace.get('peak_rss', null ) == '2.2 MB'
        trace.get('peak_vmem', null ) == '20.6 MB'

    }


}
