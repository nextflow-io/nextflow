/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.cloud.google.batch.client

import java.util.concurrent.TimeUnit

import com.google.api.core.ApiFuture
import com.google.api.core.ApiFutures
import com.google.api.gax.grpc.GrpcCallContext
import com.google.api.gax.grpc.GrpcStatusCode
import com.google.api.gax.rpc.ApiCallContext
import com.google.api.gax.rpc.Callables
import com.google.api.gax.rpc.InvalidArgumentException
import com.google.api.gax.rpc.UnaryCallable
import com.google.cloud.batch.v1.BatchServiceClient
import com.google.cloud.batch.v1.ListTasksRequest
import com.google.cloud.batch.v1.ListTasksResponse
import com.google.cloud.batch.v1.Task
import com.google.cloud.batch.v1.TaskGroupName
import com.google.cloud.batch.v1.TaskName
import com.google.cloud.batch.v1.TaskStatus
import com.google.cloud.batch.v1.stub.BatchServiceStub
import com.google.cloud.batch.v1.stub.BatchServiceStubSettings
import io.grpc.Status
import nextflow.cloud.google.GoogleOpts
import spock.lang.Specification
import spock.lang.Unroll

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

    @Unroll
    def 'should list all tasks of a job with #COUNT tasks' () {
        given:
        def project = 'project-id'
        def location = 'location-id'
        def jobId = 'job-id'
        def stub = new FakeBatchServiceStub(TaskGroupName.of(project, location, jobId, 'group0'), COUNT)
        def client = new BatchClient(projectId: project, location: location, config: new GoogleOpts([:]), batchServiceClient: BatchServiceClient.create(stub))

        when:
        def tasks = client.listTasks(jobId).toList()

        then:
        tasks.size() == COUNT
        tasks*.name.toSet().size() == COUNT
        stub.requests.size() == PAGES

        where:
        COUNT | PAGES
        1     | 1
        499   | 1
        500   | 2
        501   | 2
        1500  | 4
    }

    /**
     * Emulates the paging behaviour of the Google Batch ListTasks API: a request with no page
     * size gets the server default of 500, the page size is encoded in the returned page token,
     * and a follow-up request whose page size does not match its token is rejected
     */
    static class FakeBatchServiceStub extends BatchServiceStub {
        static final int DEFAULT_PAGE_SIZE = 500

        final List<ListTasksRequest> requests = []
        private final TaskGroupName parent
        private final int count

        FakeBatchServiceStub(TaskGroupName parent, int count) {
            this.parent = parent
            this.count = count
        }

        @Override
        UnaryCallable<ListTasksRequest, ListTasksResponse> listTasksCallable() {
            return new UnaryCallable<ListTasksRequest, ListTasksResponse>() {
                @Override
                ApiFuture<ListTasksResponse> futureCall(ListTasksRequest request, ApiCallContext context) {
                    return ApiFutures.immediateFuture(listTasks(request))
                }
            }
        }

        @Override
        UnaryCallable<ListTasksRequest, BatchServiceClient.ListTasksPagedResponse> listTasksPagedCallable() {
            return Callables
                .paged(listTasksCallable(), BatchServiceStubSettings.newBuilder().build().listTasksSettings())
                .withDefaultCallContext(GrpcCallContext.createDefault())
        }

        private ListTasksResponse listTasks(ListTasksRequest request) {
            requests.add(request)
            int offset = 0
            int pageSize = request.pageSize ?: DEFAULT_PAGE_SIZE
            if( request.pageToken ) {
                final token = request.pageToken.tokenize(':')
                offset = token[0] as int
                final tokenPageSize = token[1] as int
                if( request.pageSize != tokenPageSize )
                    throw new InvalidArgumentException("pagesize field is invalid. mismatching token page size error: request page size (${request.pageSize}) != token page size (${tokenPageSize})", null, GrpcStatusCode.of(Status.Code.INVALID_ARGUMENT), false)
            }
            final end = Math.min(offset + pageSize, count)
            final result = ListTasksResponse.newBuilder()
            for( int i = offset; i < end; i++ )
                result.addTasks(Task.newBuilder().setName("${parent}/tasks/${i}"))
            // like the real API, a full page always returns a next page token, even when no tasks remain
            if( end - offset == pageSize )
                result.setNextPageToken("${end}:${pageSize}")
            return result.build()
        }

        @Override void close() {}
        @Override void shutdown() {}
        @Override boolean isShutdown() { true }
        @Override boolean isTerminated() { true }
        @Override void shutdownNow() {}
        @Override boolean awaitTermination(long duration, TimeUnit unit) { true }
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
