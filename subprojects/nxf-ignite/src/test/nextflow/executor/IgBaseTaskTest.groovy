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

package nextflow.executor

import nextflow.processor.TaskId
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgBaseTaskTest extends Specification {

    static class DummyBaseTask extends IgBaseTask {

        DummyBaseTask(taskId) {
            super()
            def field = IgBaseTask.getDeclaredField('taskId')
            field.setAccessible(true)
            field.set(this,TaskId.of(taskId))
        }

        @Override
        protected Object execute0() {
            return null
        }

        @Override
        void cancel() {

        }
    }


    def 'should validate equals and hashCode contracts' () {

        given:
        def task_x = new DummyBaseTask(1)
        def task_y = new DummyBaseTask(1)
        def task_z = new DummyBaseTask(2)

        expect:
        task_x == task_y
        task_x != task_z

        task_x.hashCode() == task_y.hashCode()
        task_x.hashCode() != task_z.hashCode()

    }

}
