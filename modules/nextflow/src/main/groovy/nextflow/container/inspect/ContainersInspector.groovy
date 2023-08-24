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
 *
 */

package nextflow.container.inspect

import groovy.json.JsonBuilder
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.dag.DAG
import nextflow.exception.AbortOperationException
import org.codehaus.groovy.util.ListHashMap
/**
 * Preview the list of containers used by a pipeline.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ContainersInspector {

    private DAG dag

    private String format

    private boolean ignoreErrors

    ContainersInspector(DAG dag) {
        this.dag = dag
    }

    ContainersInspector withFormat(String format) {
        if( format !in ['config', 'json'] )
            throw new AbortOperationException("Invalid format for containers inspect '${format}' -- should be 'config' or 'json'")
        this.format = format
        return this
    }

    ContainersInspector withIgnoreErrors(boolean ignore) {
        this.ignoreErrors = ignore
        return this
    }

    String renderContainers() {
        log.debug "Rendering container preview"
        final containers = getContainers()
        if( format == 'config' )
            return renderConfig(containers)
        if( format == 'json' )
            return renderJson(containers)
        else
            throw new IllegalStateException("Unknown containers preview format: $format")
    }

    void printContainers() {
        final result = renderContainers()
        if( result )
            print result
    }

    protected Map<String,String> getContainers() {
        final containers = new ListHashMap<String,String>()

        for( def vertex : dag.vertices ) {
            // skip nodes that are not processes
            final process = vertex.process
            if( !process )
                continue

            try {
                // get container preview
                containers[process.name] = process.createTaskPreview().getContainer()
            }
            catch( Exception e ) {
                if( ignoreErrors )
                    log.warn "Unable to inspect container for task `$process.name` - cause: ${e.message}"
                else
                    throw e
            }
        }

        return containers
    }

    protected String renderConfig(Map<String,String> containers) {
        final result = new StringBuilder()
        for( Map.Entry<String,String> entry : containers ) {
            result.append("process { withName: '${entry.key}' { container = '${entry.value}' } }\n")
        }
        return result.toString()
    }

    protected String renderJson(Map<String,String> containers) {
        final list = containers.collect( (k, v) -> [name: k, container: v] )
        return JsonOutput.prettyPrint(new JsonBuilder(list).toString()) + '\n'
    }

}
