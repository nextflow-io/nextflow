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

package nextflow.trace
import java.nio.file.Files

import nextflow.executor.NopeTaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceFileObserverTest extends Specification {


    def testFileTrace() {

        given:
        def testFolder = Files.createTempDirectory('trace-dir')
        def file = testFolder.resolve('trace')

        // the handler
        def task = new TaskRun(id:111, name:'simple_task')
        def handler = new NopeTaskHandler(task)
        def now = System.currentTimeMillis()

        // the observer class under test
        def observer = new TraceFileObserver(file)

        when:
        observer.onFlowStart(null)
        then:
        observer.current == [:]

        when:
        observer.onProcessSubmit( handler )
        def record = observer.current.get(111)
        then:
        observer.delim == ','
        record.id == 111
        record.name == 'simple_task'
        record.submit >= now
        record.start == 0
        record.complete == 0
        observer.current.containsKey(111)

        when:
        sleep 50
        observer.onProcessStart( handler )
        then:
        record.start >= record.submit
        record.complete == 0
        observer.current.containsKey(111)

        when:
        sleep 50
        observer.onProcessComplete(handler)
        then:
        record.start >= record.submit
        record.complete >= record.start
        !(observer.current.containsKey(111))

        when:
        observer.onFlowComplete()
        def head = file.text.readLines()[0].trim()
        def parts = file.text.readLines()[1].split(',') as List
        then:
        head == observer.HEAD.join(',')
        parts[0] == '111'
        parts[1] == 'simple_task'
        TraceFileObserver.FMT.parse(parts[2]).time == record.submit
        TraceFileObserver.FMT.parse(parts[3]).time == record.start
        TraceFileObserver.FMT.parse(parts[4]).time == record.complete
        parts[5].toLong() == record.complete -record.submit
        parts[6].toLong() == record.complete -record.start


        cleanup:
        testFolder.deleteDir()

    }

}
