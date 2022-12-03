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
import com.google.cloud.logging.Severity
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.cloud.google.batch.client.BatchConfig
/**
 * Batch logging client
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class BatchLogging {
    private LoggingOptions opts
    private String projectId

    /** only for testing - do not use */
    protected BatchLogging() {

    }

    BatchLogging(BatchConfig config) {
        final creds = config.googleOpts.credentials
        this.projectId = config.googleOpts.projectId
        this.opts = LoggingOptions .newBuilder() .setCredentials(creds) .setProjectId(this.projectId) .build()
    }

    String stdout(String jobId) {
        return safeLogs(jobId,0)
    }

    String stderr(String jobId) {
        return safeLogs(jobId,1)
    }

    protected String safeLogs(String jobId, int index) {
        try {
            return fetchLogs(jobId)[index]
        }
        catch (Exception e) {
            log.warn("Cannot read logs for Batch job '$jobId' - cause: ${e.message}", e)
            return null
        }
    }

    @Memoized(maxCacheSize = 1000)
    @PackageScope List<String> fetchLogs(String uid) {
        final stdout = new StringBuilder()
        final stderr = new StringBuilder()
        try(Logging logging = opts.getService()) {
            // use logging here
            final filter = "resource.type=generic_task AND logName=\"projects/${projectId}/logs/batch_task_logs\" AND labels.job_uid=$uid"
            final entries = logging.listLogEntries(
                    Logging.EntryListOption.filter(filter),
                    Logging.EntryListOption.pageSize(1000) )

            final page = entries.getValues()
            for (LogEntry logEntry : page.iterator()) {
                parseOutput(logEntry, stdout, stderr)
            }
        }
        return [ stdout.toString(), stderr.toString() ]
    }

    protected void parseOutput(LogEntry logEntry, StringBuilder stdout, StringBuilder stderr) {
        final output = logEntry.payload.data.toString()
        if (logEntry.severity == Severity.ERROR) {
            stderr.append(output)
        } else {
            stdout.append(output)
        }
    }
}
