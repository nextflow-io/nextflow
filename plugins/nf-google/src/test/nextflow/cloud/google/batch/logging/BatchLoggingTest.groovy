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

package nextflow.cloud.google.batch.logging

import java.util.concurrent.TimeUnit

import com.google.cloud.batch.v1.Job
import com.google.cloud.batch.v1.LogsPolicy
import com.google.cloud.batch.v1.Runnable
import com.google.cloud.batch.v1.TaskGroup
import com.google.cloud.batch.v1.TaskSpec
import com.google.cloud.logging.LogEntry
import com.google.cloud.logging.Payload.StringPayload
import com.google.cloud.logging.Severity
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.cloud.google.batch.client.BatchClient
import nextflow.cloud.google.batch.client.BatchConfig
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Timeout
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BatchLoggingTest extends Specification {

    def 'should parse stdout and stderr' () {
        given:
        def OUT_ENTRY1 = LogEntry.newBuilder(StringPayload.of('No user sessions are running outdated binaries.\n')).setSeverity(Severity.INFO).build()
        def OUT_ENTRY2 = LogEntry.newBuilder(StringPayload.of('Hello world')).setSeverity(Severity.INFO).build()
        def ERR_ENTRY1 = LogEntry.newBuilder(StringPayload.of('Oops something has failed. We are sorry.\n')).setSeverity(Severity.ERROR).build()
        def ERR_ENTRY2 = LogEntry.newBuilder(StringPayload.of('blah blah')).setSeverity(Severity.ERROR).build()
        and:
        def client = new BatchLogging()

        when:
        def stdout = new StringBuilder()
        def stderr = new StringBuilder()
        and:
        client.parseOutput(OUT_ENTRY1, stdout, stderr)
        then:
        stdout.toString() == 'No user sessions are running outdated binaries.\n'
        and:
        stderr.toString() == ''

        when:
        client.parseOutput(ERR_ENTRY1, stdout, stderr)
        then:
        stderr.toString() == 'Oops something has failed. We are sorry.\n'

        when:
        client.parseOutput(ERR_ENTRY2, stdout, stderr)
        then:
        // the message is appended to the stderr because not prefix is provided
        stderr.toString() == 'Oops something has failed. We are sorry.\nblah blah'
        and:
        // no change to the stdout
        stdout.toString() == 'No user sessions are running outdated binaries.\n'

        when:
        client.parseOutput(OUT_ENTRY2, stdout, stderr)
        then:
        // the message is added to the stdout
        stdout.toString() == 'No user sessions are running outdated binaries.\nHello world'
        and:
        // no change to the stderr
        stderr.toString() == 'Oops something has failed. We are sorry.\nblah blah'

    }

    @Timeout(value = 10, unit = TimeUnit.MINUTES)
    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should fetch logs' () {
        given:
        def sess = Mock(Session) { getConfig() >> [:] }
        def config = BatchConfig.create(sess)
        and:
        def batchClient = new BatchClient(config)
        def logClient = new BatchLogging(config)

        when:
        def imageUri = 'quay.io/nextflow/bash'
        def cmd = ['/bin/bash','-c','echo "Hello world!" && echo "Oops something went wrong" >&2']
        def req = Job.newBuilder()
            .addTaskGroups(
                TaskGroup.newBuilder()
                    .setTaskSpec(
                        TaskSpec.newBuilder()
                            .addRunnables(
                                Runnable.newBuilder()
                                    .setContainer(
                                        Runnable.Container.newBuilder()
                                            .setImageUri(imageUri)
                                            .addAllCommands(cmd)
                                    )
                            )
                    )
            )
            .setLogsPolicy(
                LogsPolicy.newBuilder()
                    .setDestination(LogsPolicy.Destination.CLOUD_LOGGING)
            )
            .build()
        def jobId = 'nf-test-' + System.currentTimeMillis()
        def resp = batchClient.submitJob(jobId, req)
        def uid = resp.getUid()
        log.debug "Test job uid=$uid"
        then:
        uid
        
        when:
        def state=null
        do {
            if( batchClient.listTasks(jobId).iterator().hasNext() )
                state = batchClient.getTaskState(jobId, '0')
            else
                state = 'PENDING'
            log.debug "Test task state=$state"
            sleep 10_000
        } while( state !in ['SUCCEEDED', 'FAILED'] )
        then:
        state in ['SUCCEEDED', 'FAILED']

        when:
        def stdout = logClient.stdout(uid, '0')
        def stderr = logClient.stderr(uid, '0')
        log.debug "STDOUT: $stdout"
        log.debug "STDERR: $stderr"
        then:
        stdout.contains('Hello world!')
        stderr.contains('Oops something went wrong')
    }

}
