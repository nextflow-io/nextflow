/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */
package nextflow.cloud.google.batch.client

import com.google.cloud.batch.v1.Task
import com.google.cloud.batch.v1.TaskName
import com.google.cloud.batch.v1.TaskStatus
import spock.lang.Specification

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class BatchClientTest extends Specification{

    def 'should return task status with getTaskInArray' () {
        given:
        def project = 'project-id'
        def location = 'location-id'
        def job1 = 'job1-id'
        def task1 = 'task1-id'
        def task1Name = TaskName.of(project, location, job1, 'group0', task1).toString()
        def job2 = 'job2-id'
        def task2 = 'task2-id'
        def task2Name = TaskName.of(project, location, job2, 'group0', task2).toString()
        def job3 = 'job3-id'
        def task3 = 'task3-id'
        def task3Name = TaskName.of(project, location, job3, 'group0', task3).toString()
        def arrayTasks = new HashMap<String,TaskStatusRecord>()
        def client = Spy( new BatchClient( projectId: project, location: location, arrayTaskStatus: arrayTasks ) )

        when:
        client.listTasks(job2) >> {
            def list = new LinkedList<>()
            list.add(makeTask(task2Name, TaskStatus.State.FAILED))
            return list
        }
        client.listTasks(job3) >> {
            def list = new LinkedList<>()
            list.add(makeTask(task3Name, TaskStatus.State.SUCCEEDED))
            return list
        }
        arrayTasks.put(task1Name, makeTaskStatusRecord(TaskStatus.State.RUNNING, System.currentTimeMillis()))
        arrayTasks.put(task2Name, makeTaskStatusRecord(TaskStatus.State.PENDING, System.currentTimeMillis() - 1_001))

        then:
        // recent cached task
        client.getTaskInArrayStatus(job1, task1).state == TaskStatus.State.RUNNING
        // Outdated cached task
        client.getTaskInArrayStatus(job2, task2).state == TaskStatus.State.FAILED
        // no cached task
        client.getTaskInArrayStatus(job3, task3).state == TaskStatus.State.SUCCEEDED
    }

    TaskStatusRecord makeTaskStatusRecord(TaskStatus.State state, long timestamp) {
        return new TaskStatusRecord(TaskStatus.newBuilder().setState(state).build(), timestamp)
    }

    def makeTask(String name, TaskStatus.State state){
        Task.newBuilder().setName(name)
            .setStatus(TaskStatus.newBuilder().setState(state).build())
            .build()
    }

}
