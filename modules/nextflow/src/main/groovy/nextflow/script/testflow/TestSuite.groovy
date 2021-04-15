package nextflow.script.testflow

import java.nio.file.Path
import java.time.Duration
import java.time.LocalDateTime

/**
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class TestSuite {

    String name

    int tests
    int skipped
    int failures
    int errors

    List<TestCase> testcase

    String systemOut
    String systemErr

    Duration time
    LocalDateTime timestamp

    Path workDir

}
