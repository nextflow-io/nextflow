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

package nextflow.ga4gh.tes.server.verticle

import com.google.common.hash.HashCode
import nextflow.Session
import nextflow.ga4gh.tes.server.model.TesExecutor
import nextflow.ga4gh.tes.server.model.TesResources
import nextflow.ga4gh.tes.server.model.TesTask
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskServiceApiImplTest extends Specification {

    def 'should submit a task' () {
        given:
        def session = new Session()
        def service = new TaskServiceApiImpl(session)

        def task = new TesTask()
        def exec = new TesExecutor()
        exec.command = [ 'echo hello' ]
        exec.image = 'busybox'
        exec.workdir = '/work'
        task.id = 'task1'
        task.description = 'fake task'
        task.name = 'task-one'
        task.inputs = [ ]
        task.outputs = []
        task.resources = new TesResources()
        task.executors = [exec]

        when:
        def response = service.submitTask(task)
        sleep 1000
        session.await()
        def hash = HashCode.fromString(response.id)
        def entry = session.cache.getTaskEntry(hash, null)
        //println entry.trace

        then:
        entry.trace.status == 'COMPLETED'
        entry.trace.name == 'task-one'
        entry.trace.exit == 0
        entry.trace.script == 'echo hello'
        entry.trace.container == 'busybox'
        //entry.trace.process == 'todo'

    }
}
