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
