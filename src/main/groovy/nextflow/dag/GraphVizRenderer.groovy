/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.dag
import java.nio.file.Files
import java.nio.file.Path

import groovy.util.logging.Slf4j
/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
@Slf4j
class GraphvizRenderer implements DagRenderer {

    private String format

    private String name

    GraphvizRenderer(String name, String format) {
        this.name = name
        this.format = format
    }

    /**
     * Render the DAG using Graphviz to the specified
     * file in a format specified in the constructor.
     * See http://www.graphviz.org for more info.
     */
    @Override
    void renderDocument(DAG dag, Path file) {
        def temp = Files.createTempFile('nxf-','.dot')
        // save the DAG as `dot` to a temp file
        temp.text = new DotRenderer(name).renderNetwork(dag)

        final cmd = "command -v dot &>/dev/null || exit 128 && dot -T${format} ${temp} > ${file}"
        final exitStatus = ["bash","-c", cmd].execute().waitFor()
        if( exitStatus == 128 ) {
            file = file.resolveSibling( "${file.baseName}.dot" )
            temp.moveTo(file)
            log.warn "To render the execution DAG in the required format it is required to install Graphviz -- See http://www.graphviz.org for more info."
        }
        else if( exitStatus>0 ) {
            log.debug("Graphviz error -- command `$cmd` -- exist status: $exitStatus")
            log.warn "Failed to render DAG file: $file"
        }
    }
}
