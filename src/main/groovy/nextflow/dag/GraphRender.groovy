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
@CompileStatic
class GraphRender implements TraceObserver {

    static final String DEF_FILE_NAME = 'dag.dot'

    private Path file

    private DAG dag

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

        final format = file.getExtension() ?: 'dot'

        if( format == 'html' ) {

           file.text = CytoscapeJsRenderer.render(dag)

        } else {

            def temp = Files.createTempFile('nxf-','.dot')
            // save the DAG as `dot` to a temp file
            temp.text = DotRenderer.render(dag)

            if( format == 'dot' ) {
                temp.moveTo(file)
            }
            else {
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
