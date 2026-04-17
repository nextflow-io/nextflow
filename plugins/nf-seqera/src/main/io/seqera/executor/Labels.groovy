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

package io.seqera.executor

import com.google.common.hash.Hashing

import groovy.transform.CompileStatic
import nextflow.NextflowMeta
import nextflow.script.WorkflowMetadata

/**
 * Helper class to manage run labels.
 *
 * Builds the labels map from workflow metadata ({@code nextflow.io/*}),
 * scheduler metadata ({@code seqera:sched:*}), and user-configured labels.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class Labels {

    private final Map<String,String> entries = new LinkedHashMap<>(20)

    /**
     * Add {@code nextflow.io/*} labels from workflow metadata
     */
    Labels withWorkflowMetadata(WorkflowMetadata workflow) {
        if( workflow.projectName )
            entries.put('nextflow.io/projectName', workflow.projectName)
        if( workflow.userName )
            entries.put('nextflow.io/userName', workflow.userName)
        if( workflow.runName )
            entries.put('nextflow.io/runName', workflow.runName)
        if( workflow.sessionId )
            entries.put('nextflow.io/sessionId', workflow.sessionId.toString())
        entries.put('nextflow.io/resume', String.valueOf(workflow.resume))
        if( workflow.revision )
            entries.put('nextflow.io/revision', workflow.revision)
        if( workflow.commitId )
            entries.put('nextflow.io/commitId', workflow.commitId)
        if( workflow.repository )
            entries.put('nextflow.io/repository', workflow.repository)
        if( workflow.manifest?.name )
            entries.put('nextflow.io/manifestName', workflow.manifest.name)
        if( NextflowMeta.instance.version )
            entries.put('nextflow.io/runtimeVersion', NextflowMeta.instance.version.toString())
        if( workflow.platform?.workflowId )
            entries.put('seqera.io/platform/workflowId', workflow.platform.workflowId)
        return this
    }

    /**
     * Add {@code seqera:sched:*} scheduler labels
     */
    Labels withSchedRunId(String runId) {
        if( runId )
            entries.put('seqera:sched:runId', runId)
        return this
    }

    Labels withSchedClusterId(String clusterId) {
        if( clusterId )
            entries.put('seqera:sched:clusterId', clusterId)
        return this
    }

    /**
     * Add config-level {@code process.resourceLabels}. Values are coerced to
     * string via {@link String#valueOf} to satisfy the scheduler API typing.
     */
    Labels withProcessResourceLabels(Map<String,?> map) {
        if( !map ) return this
        for( Map.Entry<String,?> entry : map.entrySet() )
            entries.put(entry.key.toString(), String.valueOf(entry.value))
        return this
    }

    /**
     * @return all labels as an unmodifiable map
     */
    Map<String,String> getEntries() {
        return Collections.unmodifiableMap(entries)
    }

    /**
     * Compute a run identifier as SipHash of sessionId + runName
     */
    protected static String runId(String sessionId, String runName) {
        return Hashing
                .sipHash24()
                .newHasher()
                .putUnencodedChars(sessionId)
                .putUnencodedChars(runName)
                .hash()
                .toString()
    }

    /**
     * Coerce arbitrary map values to strings via {@link String#valueOf}.
     * Returns an empty map for null/empty input.
     */
    static Map<String,String> toStringMap(Map<String,?> map) {
        if( !map ) return Collections.<String,String>emptyMap()
        final result = new LinkedHashMap<String,String>(map.size())
        for( Map.Entry<String,?> entry : map.entrySet() )
            result.put(entry.key.toString(), String.valueOf(entry.value))
        return result
    }

    /**
     * Return the entries of {@code task} that are missing from {@code run}
     * or have a different value. Returns {@code null} if the resulting
     * map would be empty (so callers can omit the field).
     */
    static Map<String,String> delta(Map<String,String> task, Map<String,String> run) {
        if( !task ) return null
        final result = new LinkedHashMap<String,String>()
        for( Map.Entry<String,String> entry : task.entrySet() ) {
            final k = entry.key
            final v = entry.value
            if( run == null || !run.containsKey(k) || run.get(k) != v )
                result.put(k, v)
        }
        return result.isEmpty() ? null : result
    }
}
