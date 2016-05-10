

package nextflow.dag

import groovy.util.logging.Slf4j
import java.nio.file.Path

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Mike Smoot <mes@aescon.com>
 */
@Slf4j
class ErrorRenderer {

    void renderDocument(DAG dag, Path file) {
        log.error "Failed to render DAG file: $file. Unrecognized file format."
    }
}
