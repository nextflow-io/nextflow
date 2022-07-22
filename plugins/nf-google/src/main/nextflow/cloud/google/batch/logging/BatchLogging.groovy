/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import com.google.cloud.logging.LogEntry
import com.google.cloud.logging.Logging
import com.google.cloud.logging.LoggingOptions
import groovy.transform.Memoized
import groovy.transform.PackageScope
import nextflow.cloud.google.batch.client.BatchConfig
/**
 * Batch logging client
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchLogging {

    private static final String STDOUT = 'STDOUT:  '
    private static final String STDERR = 'STDERR:  '

    private String mode = STDOUT
    private LoggingOptions opts
    private String projectId

    /** only for testing - do not use */
    protected BatchLogging() {

    }

    BatchLogging(BatchConfig config) {
        final creds = config.googleOpts.credentials
        this.opts = LoggingOptions .newBuilder() .setCredentials(creds) .build()
        this.projectId = config.googleOpts.projectId
    }

    String stdout(String jobId) {
        return fetchLogs(jobId)[0]
    }

    String stderr(String jobId) {
        return fetchLogs(jobId)[1]
    }

    @PackageScope String currentMode() { mode }

    @Memoized(maxCacheSize = 1000)
    @PackageScope List<String> fetchLogs(String uid) {
        try(Logging logging = opts.getService()) {
            // use logging here
            final filter = "resource.type=generic_task AND logName=projects/${projectId}/logs/batch_task_logs AND labels.job_uid=$uid"
            final entries = logging.listLogEntries(
                    Logging.EntryListOption.filter(filter),
                    Logging.EntryListOption.pageSize(1000) )

            final stdout = new StringBuilder()
            final stderr = new StringBuilder()
            final page = entries.getValues()
            for (LogEntry logEntry : page.iterator()) {
                final output = logEntry.payload.data.toString()
                parseOutput(output, stdout, stderr)
            }

            return [ stdout.toString(), stderr.toString() ]
        }
    }

    protected void parseOutput(String output, StringBuilder stdout, StringBuilder stderr) {
        // check stderr
        def p = output.indexOf(STDERR)
        if( p>=0 ) {
            mode = STDERR
            output = output.substring(p+STDERR.size())
        }
        else if( (p = output.indexOf(STDOUT))>=0 )  {
            mode = STDOUT
            output = output.substring(p+STDOUT.size())
        }
        // now append the result
        if( mode==STDOUT )
            stdout.append(output)
        else
            stderr.append(output)
    }

}
