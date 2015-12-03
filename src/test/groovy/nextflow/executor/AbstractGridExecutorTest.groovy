/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractGridExecutorTest extends Specification {

    def 'should remove invalid chars from name' () {

        given:
        def task = new TaskRun(name: 'task 90 (foo:bar)')
        def exec = [:] as AbstractGridExecutor

        expect:
        exec.getJobNameFor(task) == 'nf-task_90_(foo_bar)'

    }

    def 'should return the kill list' () {

        given:
        def exec = [getKillCommand: { 'qdel' }] as AbstractGridExecutor

        expect:
        exec.killTaskCommand('10') == ['qdel', '10']
        exec.killTaskCommand([11,12]) == ['qdel', '11', '12']
        exec.killTaskCommand([100,200,300]) == ['qdel', '100', '200', '300']

    }
}
