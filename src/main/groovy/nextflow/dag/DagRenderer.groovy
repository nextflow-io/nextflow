

package nextflow.dag
import java.nio.file.Path

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
interface DagRenderer {

    /**
     * Render the dag to the specified file.
     */
    void renderDocument(DAG dag, Path file);
}
