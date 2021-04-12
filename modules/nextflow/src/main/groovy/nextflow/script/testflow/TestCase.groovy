package nextflow.script.testflow

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

}
