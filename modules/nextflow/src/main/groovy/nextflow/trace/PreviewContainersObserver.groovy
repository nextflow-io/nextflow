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

import java.nio.file.Path

import groovy.json.JsonBuilder
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
/**
 * Render the preview report when running a pipeline
 * in preview mode.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PreviewContainersObserver implements TraceObserver {

    static final String DEF_FILE_NAME = "preview-${TraceHelper.launchTimestampFmt()}.json"

    private DAG dag

    private Path file

    private String format

    PreviewContainersObserver(Path file ) {
        this.file = file
        this.format = file.getExtension().toLowerCase() ?: 'json'

        if( format !in ['config', 'json'] )
            throw new IllegalArgumentException("Invalid format for container preview: '${format}' -- should be *.json or *.config")
    }

    @Override
    void onFlowCreate(Session session) {
        this.dag = session.dag
    }

    @Override
    void onFlowBegin() {
        log.debug "Rendering container preview"
        try {
            final containers = getContainers()
            if( format == 'config' )
                renderConfig(containers)
            else if( format == 'json' )
                renderJson(containers)
        }
        catch( Exception e ) {
            log.warn "Failed to preview containers -- see the log file for details", e
        }
    }

    private Map<String,String> getContainers() {
        final containers = [:]

        for( def vertex : dag.vertices ) {
            // skip nodes that are not processes
            final process = vertex.process
            if( !process )
                continue

            // get container preview
            try {
                containers[process.name] = process.getPreviewTask().getContainer()
            }
            catch( Exception e ) {
                log.warn1 "Unable to preview container for process `${process.name}`: ${e}"
            }
        }

        return containers
    }

    private void renderConfig(Map<String,String> containers) {
        final result = new StringBuilder()
        for( def entry : containers ) {
            result.append("process { withName: '${entry.key}' { container = '${entry.value}' } }\n")
        }

        file.text = result.toString()
    }

    private void renderJson(Map<String,String> containers) {
        file.text = JsonOutput.prettyPrint(new JsonBuilder(containers).toString())
    }

}
