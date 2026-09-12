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

package nextflow.platform

import groovy.transform.CompileStatic
import nextflow.NextflowMeta
import nextflow.script.WorkflowMetadata

/**
 * Derive resource labels from the workflow metadata.
 *
 * The labels are selected by short name (e.g. {@code runName}) via the
 * {@code tower.autoLabels} option — or the deprecated {@code seqera.executor.autoLabels} —
 * and mapped to the canonical {@code nextflow.io/*} and {@code seqera.io/platform/*} keys.
 *
 * Entries whose source metadata is absent are omitted, therefore a run without
 * Seqera Platform still gets the labels the runtime knows on its own.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AutoLabels {

    /**
     * The workflow metadata short names that can be turned into resource labels
     */
    static final Set<String> VALID_NAMES = Collections.unmodifiableSet(new LinkedHashSet<>([
        'projectName', 'userName', 'runName', 'sessionId', 'resume',
        'revision', 'commitId', 'repository', 'manifestName',
        'runtimeVersion', 'workflowId', 'workspaceId', 'computeEnvId'
    ]))

    /**
     * Parse the value of an auto-labels config option.
     *
     * @param value
     *      The config value: {@code null} or {@code false} (none), {@code true} (all names),
     *      a list of short names, or a comma-separated string of short names
     * @param optionName
     *      The name of the config option being parsed, used to make the error message actionable
     * @return
     *      The set of selected short names, empty when the feature is disabled
     */
    static Set<String> parse(Object value, String optionName='autoLabels') {
        if( value == null || value == false )
            return Collections.<String>emptySet()
        if( value == true )
            return VALID_NAMES
        List<String> raw
        if( value instanceof CharSequence )
            raw = value.toString().tokenize(',').collect { String s -> s.trim() }.findAll { String s -> s }
        else if( value instanceof List )
            raw = ((List) value).collect { it?.toString()?.trim() }.findAll { String s -> s } as List<String>
        else
            throw new IllegalArgumentException("Invalid '${optionName}' value '${value}' - expected true, false, a list, or a comma-separated string")
        final invalid = raw.findAll { String s -> !(s in VALID_NAMES) }
        if( invalid )
            throw new IllegalArgumentException("Invalid '${optionName}' name(s) ${invalid} - valid names are: ${VALID_NAMES.join(', ')}")
        return Collections.unmodifiableSet(new LinkedHashSet<>(raw))
    }

    /**
     * Map the workflow metadata to the corresponding resource labels.
     *
     * @param workflow
     *      The current {@link WorkflowMetadata} object
     * @param include
     *      The set of short names to be included, as returned by {@link #parse}
     * @return
     *      The labels map, holding only the entries whose metadata value is available
     */
    static Map<String,String> labelsFor(WorkflowMetadata workflow, Set<String> include) {
        if( !workflow || !include )
            return Collections.<String,String>emptyMap()
        final entries = new LinkedHashMap<String,String>(20)
        if( include.contains('projectName') && workflow.projectName )
            entries.put('nextflow.io/projectName', workflow.projectName)
        // the Platform user that submitted the run, falling back to the OS user running the launcher
        final userName = workflow.platform?.user?.userName ?: workflow.userName
        if( include.contains('userName') && userName )
            entries.put('nextflow.io/userName', userName)
        if( include.contains('runName') && workflow.runName )
            entries.put('nextflow.io/runName', workflow.runName)
        if( include.contains('sessionId') && workflow.sessionId )
            entries.put('nextflow.io/sessionId', workflow.sessionId.toString())
        if( include.contains('resume') )
            entries.put('nextflow.io/resume', String.valueOf(workflow.resume))
        if( include.contains('revision') && workflow.revision )
            entries.put('nextflow.io/revision', workflow.revision)
        if( include.contains('commitId') && workflow.commitId )
            entries.put('nextflow.io/commitId', workflow.commitId)
        if( include.contains('repository') && workflow.repository )
            entries.put('nextflow.io/repository', workflow.repository)
        if( include.contains('manifestName') && workflow.manifest?.name )
            entries.put('nextflow.io/manifestName', workflow.manifest.name)
        if( include.contains('runtimeVersion') && NextflowMeta.instance.version )
            entries.put('nextflow.io/runtimeVersion', NextflowMeta.instance.version.toString())
        if( include.contains('workflowId') && workflow.platform?.workflowId )
            entries.put('seqera.io/platform/workflowId', workflow.platform.workflowId)
        if( include.contains('workspaceId') && workflow.platform?.workspace?.id )
            entries.put('seqera.io/platform/workspaceId', workflow.platform.workspace.id)
        if( include.contains('computeEnvId') && workflow.platform?.computeEnv?.id )
            entries.put('seqera.io/platform/computeEnvId', workflow.platform.computeEnv.id)
        return Collections.unmodifiableMap(entries)
    }

}
