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
import com.google.api.gax.rpc.AlreadyExistsException
import com.google.api.gax.rpc.ApiCallContext
import com.google.api.gax.rpc.Callables
import com.google.api.gax.rpc.DeadlineExceededException
import com.google.api.gax.rpc.InvalidArgumentException
import com.google.api.gax.rpc.StatusCode
import com.google.api.gax.rpc.UnaryCallable
import com.google.cloud.batch.v1.BatchServiceClient
import com.google.cloud.batch.v1.CreateJobRequest
import com.google.cloud.batch.v1.GetJobRequest
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.JobName
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

    /**
     * A real {@link BatchServiceClient} over a stub transport: its RPC methods are final, so
     * they cannot be mocked, but the stub it delegates to can
     */
    private BatchServiceClient createService(List<Closure> createJob, Map<String,Job> jobs, List<String> getJobCalls) {
        def stub = new BatchServiceStub() {
            @Override
            UnaryCallable<CreateJobRequest, Job> createJobCallable() {
                return new UnaryCallable<CreateJobRequest, Job>() {
                    @Override
                    ApiFuture<Job> futureCall(CreateJobRequest request, ApiCallContext context) {
                        try {
                            return ApiFutures.immediateFuture((Job) createJob.remove(0).call(request))
                        }
                        catch( Throwable t ) {
                            return ApiFutures.<Job>immediateFailedFuture(t)
                        }
                    }
                }
            }
            @Override
            UnaryCallable<GetJobRequest, Job> getJobCallable() {
                return new UnaryCallable<GetJobRequest, Job>() {
                    @Override
                    ApiFuture<Job> futureCall(GetJobRequest request, ApiCallContext context) {
                        getJobCalls << request.getName()
                        return ApiFutures.immediateFuture(jobs.get(request.getName()))
                    }
                }
            }
            @Override void close() {}
            @Override void shutdown() {}
            @Override boolean isShutdown() { false }
            @Override boolean isTerminated() { false }
            @Override void shutdownNow() {}
            @Override boolean awaitTermination(long duration, TimeUnit unit) { true }
        }
        return BatchServiceClient.create(stub)
    }

    private BatchClient createClient(BatchServiceClient service) {
        def client = new BatchClient()
        client.projectId = 'project-id'
        client.location = 'location-id'
        client.config = new GoogleOpts([project: 'project-id', location: 'location-id', batch: [retryPolicy: [delay: '1ms', maxDelay: '10ms', maxAttempts: 3]]])
        client.batchServiceClient = service
        return client
    }

    def 'should reuse the job created by a previous submit attempt' () {
        given:
        def name = JobName.of('project-id', 'location-id', 'job-1').toString()
        def created = Job.newBuilder().setName(name).build()
        def getJobCalls = []
        def service = createService([
                // the first attempt creates the job server-side but the response is lost
                { throw new DeadlineExceededException(new RuntimeException('deadline'), Stub(StatusCode), true) },
                // the retry finds the job it submitted itself
                { throw new AlreadyExistsException(new RuntimeException('exists'), Stub(StatusCode), false) } ],
                [(name): created], getJobCalls)
        def client = createClient(service)

        when:
        def result = client.submitJob('job-1', Job.newBuilder().build())

        then:
        result == created
        getJobCalls == [name]
    }

    def 'should fail when the job already exists on the first submit attempt' () {
        given:
        def getJobCalls = []
        def service = createService([
                { throw new AlreadyExistsException(new RuntimeException('exists'), Stub(StatusCode), false) } ],
                [:], getJobCalls)
        def client = createClient(service)

        when:
        client.submitJob('job-1', Job.newBuilder().build())

        then:
        thrown(AlreadyExistsException)
        getJobCalls == []
    }

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
