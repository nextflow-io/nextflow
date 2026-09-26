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

import com.google.api.gax.core.NoCredentialsProvider
import com.google.api.gax.grpc.GrpcTransportChannel
import com.google.api.gax.rpc.AlreadyExistsException
import com.google.api.gax.rpc.FixedTransportChannelProvider
import com.google.api.gax.rpc.NotFoundException
import com.google.cloud.batch.v1.BatchServiceClient
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.cloud.batch.v1.CreateJobRequest
import com.google.cloud.batch.v1.GetJobRequest
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.JobName
import io.grpc.CallOptions
import io.grpc.ClientCall
import io.grpc.ManagedChannel
import io.grpc.Metadata
import io.grpc.MethodDescriptor
import io.grpc.Status
import nextflow.cloud.google.GoogleOpts
import spock.lang.Specification

/**
 * Checks that a `createJob` request whose response is lost does not abort the task submission.
 *
 * `BatchServiceClient.createJob` and `getJob` are final and cannot be stubbed, so the real
 * client is driven against a fake channel that answers calls with a scripted sequence of
 * results. This exercises the actual API client and the failsafe retry policy around it.
 *
 * see https://github.com/nextflow-io/nextflow/issues/6916
 */
class BatchClientSubmitJobTest extends Specification {

    static final String PROJECT = 'proj-1'
    static final String LOCATION = 'europe-west1'
    static final String JOB_ID = 'nf-06071afa-1773255240303'
    static final String JOB_NAME = JobName.of(PROJECT, LOCATION, JOB_ID).toString()

    /**
     * A fake Batch endpoint replaying a scripted sequence of results per RPC. Each queued entry
     * is either a {@link Status} to fail the call with, or a {@link Job} to return.
     */
    static class FakeBatchChannel extends ManagedChannel {

        final List<String> calls = Collections.synchronizedList(new ArrayList<String>())
        final Queue<Object> createResults = new LinkedList<>()
        final Queue<Object> getResults = new LinkedList<>()

        @Override
        String authority() { 'batch.googleapis.test' }

        @Override
        ClientCall newCall(MethodDescriptor method, CallOptions options) {
            return new ClientCall() {
                private ClientCall.Listener listener
                private Object result

                @Override
                void start(ClientCall.Listener responseListener, Metadata headers) {
                    this.listener = responseListener
                }

                @Override
                void request(int numMessages) { }

                @Override
                void cancel(String message, Throwable cause) { }

                @Override
                void sendMessage(Object message) {
                    // record the call and pick the next scripted result for this RPC
                    if( message instanceof CreateJobRequest ) {
                        calls << "createJob:${message.getJobId()}".toString()
                        result = createResults.poll()
                    }
                    else if( message instanceof GetJobRequest ) {
                        calls << "getJob:${message.getName()}".toString()
                        result = getResults.poll()
                    }
                    else {
                        throw new IllegalStateException("Unexpected request: ${message?.getClass()}")
                    }
                    if( result == null )
                        throw new IllegalStateException("No scripted result left for ${method.getFullMethodName()}")
                }

                @Override
                void halfClose() {
                    if( result instanceof Status ) {
                        listener.onClose(result as Status, new Metadata())
                    }
                    else {
                        listener.onMessage(result)
                        listener.onClose(Status.OK, new Metadata())
                    }
                }
            }
        }

        @Override ManagedChannel shutdown() { return this }
        @Override ManagedChannel shutdownNow() { return this }
        @Override boolean isShutdown() { return false }
        @Override boolean isTerminated() { return false }
        @Override boolean awaitTermination(long timeout, TimeUnit unit) { return true }
    }

    FakeBatchChannel fake
    BatchServiceClient service

    def setup() {
        fake = new FakeBatchChannel()
        final settings = BatchServiceSettings.newBuilder()
            .setTransportChannelProvider(FixedTransportChannelProvider.create(GrpcTransportChannel.create(fake)))
            .setCredentialsProvider(NoCredentialsProvider.create())
            .build()
        service = BatchServiceClient.create(settings)
    }

    def cleanup() {
        service?.close()
    }

    BatchClient newClient() {
        new BatchClient(
            projectId: PROJECT,
            location: LOCATION,
            batchServiceClient: service,
            config: new GoogleOpts([batch: [retryPolicy: [delay: '1ms', maxDelay: '10ms']]]) )
    }

    static Job aJob(String uid) {
        Job.newBuilder().setName(JOB_NAME).setUid(uid).build()
    }

    def 'should submit a job' () {
        given:
        fake.createResults << aJob('job-uid-0')

        when:
        def result = newClient().submitJob(JOB_ID, Job.newBuilder().build())

        then:
        result.getUid() == 'job-uid-0'
        and:
        fake.calls == ["createJob:$JOB_ID"]
    }

    def 'should adopt the job created by a request whose response was lost' () {
        given: 'the create reaches Batch but its response is lost, so the retry finds the job it created'
        fake.createResults << Status.UNAVAILABLE.withDescription('502:Bad Gateway')
        fake.createResults << Status.ALREADY_EXISTS.withDescription("Resource '$JOB_NAME' already exists")
        fake.getResults << aJob('job-uid-1')

        when:
        def result = newClient().submitJob(JOB_ID, Job.newBuilder().build())

        then: 'the submission returns the pre-existing job instead of aborting the task'
        result.getUid() == 'job-uid-1'
        and:
        fake.calls == ["createJob:$JOB_ID", "createJob:$JOB_ID", "getJob:$JOB_NAME"]
    }

    def 'should recover the same way from a client-side deadline' () {
        given:
        fake.createResults << Status.DEADLINE_EXCEEDED.withDescription('deadline exceeded after 59.999709386s')
        fake.createResults << Status.ALREADY_EXISTS.withDescription("Resource '$JOB_NAME' already exists")
        fake.getResults << aJob('job-uid-3')

        when:
        def result = newClient().submitJob(JOB_ID, Job.newBuilder().build())

        then:
        result.getUid() == 'job-uid-3'
    }

    def 'should retry the create when the job cannot be read back' () {
        given: 'ALREADY_EXISTS is reported but the job is not there'
        fake.createResults << Status.UNAVAILABLE.withDescription('502:Bad Gateway')
        fake.createResults << Status.ALREADY_EXISTS.withDescription("Resource '$JOB_NAME' already exists")
        fake.getResults << Status.NOT_FOUND.withDescription("Resource '$JOB_NAME' was not found")
        and: 'the next create attempt succeeds'
        fake.createResults << aJob('job-uid-2')

        when:
        def result = newClient().submitJob(JOB_ID, Job.newBuilder().build())

        then: 'the lookup failure retries the create instead of aborting the task'
        result.getUid() == 'job-uid-2'
        and:
        fake.calls == ["createJob:$JOB_ID", "createJob:$JOB_ID", "getJob:$JOB_NAME", "createJob:$JOB_ID"]
    }

    def 'should give up when the job can neither be read back nor recreated' () {
        given: 'every create is reported ALREADY_EXISTS and every lookup NOT_FOUND'
        5.times {
            fake.createResults << Status.ALREADY_EXISTS.withDescription("Resource '$JOB_NAME' already exists")
            fake.getResults << Status.NOT_FOUND.withDescription("Resource '$JOB_NAME' was not found")
        }

        when:
        newClient().submitJob(JOB_ID, Job.newBuilder().build())

        then: 'the retry policy is bounded and the error surfaces'
        thrown(NotFoundException)
        and:
        fake.calls.count { it.startsWith('createJob') } == 5
    }

    def 'should surface ALREADY_EXISTS from a generic API call' () {
        given: 'the recovery is specific to submitJob -- other calls keep the plain behaviour'
        fake.getResults << Status.ALREADY_EXISTS.withDescription('unexpected')

        when:
        newClient().describeJob(JOB_ID)

        then:
        thrown(AlreadyExistsException)
    }
}
