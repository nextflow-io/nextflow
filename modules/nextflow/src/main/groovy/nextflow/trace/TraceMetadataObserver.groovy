/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.trace

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.agent.Agent
import nextflow.ISession
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.script.WorkflowMetadata

import java.nio.file.Files
import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap

@Slf4j
@CompileStatic
class TraceMetadataObserver implements TraceObserver {

    Session session

    TraceMetadataObserver(Session session) {
        this.session = session
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        final taskId = handler.task.id
        if( !trace ) {
            log.debug "[WARN] Unable to find record for task run with id: ${taskId}"
            return
        }

        session.workflowMetadata.traces.add(trace)
    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {
        onProcessComplete(handler, trace)
    }

    @Override
    boolean enableMetrics() {
        return true
    }
}
