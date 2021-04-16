package nextflow.script.testflow

import java.nio.file.Path
import java.time.Duration

/**
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class TestCase {

    String name
    String className
    TestFailure failure
    Duration time
    Path workDir

    boolean isFailed() {
        failure != null
    }

}
