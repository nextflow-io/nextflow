/*
 * Copyright 2022, Google Inc.
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

import com.google.api.gax.core.CredentialsProvider
import com.google.auth.Credentials
import com.google.cloud.batch.v1.BatchServiceClient
import com.google.cloud.batch.v1.BatchServiceSettings
import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.JobName
import com.google.cloud.batch.v1.JobStatus
import com.google.cloud.batch.v1.LocationName
import com.google.cloud.batch.v1.TaskGroupName
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 * Implements Google Batch HTTP client
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class BatchClient {

    protected String projectId
    protected String location
    protected BatchServiceClient batchServiceClient

    BatchClient(BatchConfig config) {
        this.projectId = config.googleOpts.projectId
        this.location = config.googleOpts.location
        this.batchServiceClient = createBatchService(config)
    }

    /** Only for testing - do not use */
    protected BatchClient() {}

    protected CredentialsProvider createCredentialsProvider(BatchConfig config) {
        if( !config.getCredentials() )
            return null
        return new CredentialsProvider() {
            @Override
            Credentials getCredentials() throws IOException {
                return config.getCredentials()
            }
        }
    }

    protected BatchServiceClient createBatchService(BatchConfig config) {
        final provider = createCredentialsProvider(config)
        if( provider ) {
            log.debug "[GOOGLE BATCH] Creating service client with config credentials"
            final settings = BatchServiceSettings.newBuilder().setCredentialsProvider(provider).build()
            return BatchServiceClient.create(settings)
        }
        else {
            log.debug "[GOOGLE BATCH] Creating service client with default settings"
            return BatchServiceClient.create()
        }
    }

    Job submitJob(String jobId, Job job) {
        final parent = LocationName.of(projectId, location)

        return batchServiceClient.createJob(parent, job, jobId)
    }

    Job describeJob(String jobId) {
        final name = JobName.of(projectId, location, jobId)

        return batchServiceClient.getJob(name)
    }

    Iterable<?> listTasks(String jobId) {
        final parent = TaskGroupName.of(projectId, location, jobId, 'group0')

        return batchServiceClient.listTasks(parent).iterateAll()
    }

    void deleteJob(String jobId) {
        final name = JobName.of(projectId, location, jobId).toString()

        batchServiceClient.deleteJobAsync(name)
    }

    JobStatus getJobStatus(String jobId) {
        final job = describeJob(jobId)
        return job.getStatus()
    }

    String getJobState(String jobId) {
        final status = getJobStatus(jobId)
        return status ? status.getState().toString() : null
    }

    void shutdown() {
        batchServiceClient.close()
    }
}
