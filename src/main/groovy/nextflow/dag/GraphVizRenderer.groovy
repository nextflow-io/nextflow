/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
