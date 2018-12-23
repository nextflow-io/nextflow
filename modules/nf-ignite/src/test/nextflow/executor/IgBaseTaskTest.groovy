/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
