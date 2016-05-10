
package nextflow.dag

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.nio.file.Path
import java.nio.file.Files

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
@Slf4j
class GraphVizRenderer implements DagRenderer {

    private final String format

    GraphVizRenderer(String format) {
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
        temp.text = DotRenderer.renderNetwork(dag)

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
