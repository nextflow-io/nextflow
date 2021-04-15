package nextflow.script.testflow

import java.time.Duration

/**
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class TestRun {

    private List<TestSuite> testSuites = new ArrayList<>()

    void addSuite(TestSuite suite) {
        testSuites.add(suite)
    }

    List<TestSuite> getSuites() {
        testSuites
    }

    boolean isEmpty() {
        testSuites.isEmpty()
    }

    int getCompleted() {
        if (isEmpty()) {
            return 0
        }
        return testSuites.collect{ it.tests }.sum() as int
    }

    int getSkipped() {
        if (isEmpty()) {
            return 0
        }
        return testSuites.collect{ it.skipped }.sum() as int
    }

    int getFailures() {
        if (isEmpty()) {
            return 0
        }
        return testSuites.collect{ it.failures }.sum() as int
    }

    int getErrors() {
        if (isEmpty())
            return 0
        return testSuites.collect{ it.errors }.sum() as int
    }

    int getFailed() {
        failures + errors
    }

    Duration getDuration() {
        testSuites.collect { it.time }.sum() as Duration
    }

    int getSuccessRate() {
        ((completed - failed) * 100 / completed) as int
    }

}
