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

package nextflow.trace
import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.CytoscapeHtmlRenderer
import nextflow.dag.DAG
import nextflow.dag.DagRenderer
import nextflow.dag.DotRenderer
import nextflow.dag.GraphvizRenderer
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
/**
 * Render the DAG document on pipeline completion using the
 * format specified by the user
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class GraphObserver implements TraceObserver {

    static public final String DEF_FILE_NAME = 'dag.dot'

    private Path file

    private DAG dag

    private String name

    private String format

    String getFormat() { format }

    String getName() { name }

    GraphObserver( Path file ) {
        assert file
        file.rollFile()
        this.file = file
        this.name = file.baseName
        this.format = file.getExtension().toLowerCase() ?: 'dot'
    }

    @Override
    void onFlowStart(Session session) {
        this.dag = session.dag
    }

    @Override
    void onFlowComplete() {
        // -- normalise the DAG
        dag.normalize()
        // -- render it to a file
        createRender().renderDocument(dag,file)
    }

    @PackageScope
    DagRenderer createRender() {
        if( format == 'dot' )
            new DotRenderer(name)

        else if( format == 'html' )
            new CytoscapeHtmlRenderer()

        else
            new GraphvizRenderer(name, format)
    }


    @Override
    void onProcessCreate(TaskProcessor process) {

    }


    @Override
    void onProcessSubmit(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    void onProcessStart(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    void onProcessCached(TaskHandler handler, TraceRecord trace) {

    }

    @Override
    boolean enableMetrics() {
        return false
    }
}
