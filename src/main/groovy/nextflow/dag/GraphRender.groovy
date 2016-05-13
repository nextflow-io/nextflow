/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.trace.TraceObserver

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
//@CompileStatic
class GraphRender implements TraceObserver {

    static final String DEF_FILE_NAME = 'dag.dot'

    private Path file

    private DAG dag

    private Map<String, DagRenderer> renderers = [ 'png' : new GraphVizRenderer('png'),
                                                   'pdf' : new GraphVizRenderer('pdf'),
                                                   'svg' : new GraphVizRenderer('svg'),
                                                   'dot' : new DotRenderer(),
                                                   'htm' : new DagreD3Renderer(),
                                                   'html': new CytoscapeHtmlRenderer(),
                                                   'json': new CytoscapeJsRenderer() ]


    GraphRender( Path file ) {
        this.file = file
    }

    @Override
    void onFlowStart(Session session) {
        this.dag = session.dag
    }

    @Override
    void onFlowComplete() {
        dag.normalize()
        final format = (file.getExtension() ?: 'dot').toLowerCase()
        renderers.get(format, new ErrorRenderer()).renderDocument(dag, file)
    }


    @Override
    void onProcessCreate(TaskProcessor process) {

    }


    @Override
    void onProcessSubmit(TaskHandler handler) {

    }

    @Override
    void onProcessStart(TaskHandler handler) {

    }

    @Override
    void onProcessComplete(TaskHandler handler) {

    }

}
