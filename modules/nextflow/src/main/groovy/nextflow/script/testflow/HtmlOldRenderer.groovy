package nextflow.script.testflow

import com.kncept.junit.reporter.TestReportProcessor
import com.kncept.junit.reporter.TestRunResults
import com.kncept.junit.reporter.domain.CssRagStatus
import com.kncept.junit.reporter.exception.TestReporterFailure
import com.kncept.junit.reporter.logger.Log
import com.kncept.junit.reporter.logger.LogFactory

import java.nio.file.Path

/**
 * Creates an HTML test report from the XUnit XML files at the test-results folder
 *
 * @author Jordi Deu-Pons <jordi@jordeu.net>
 */
class HtmlOldRenderer {

    /**
     * Write test HTML report
     *
     * @param testDir Folder that contains all the testsuite XUnit XML files
     * @param reportDir Output folder to write the HTML
     */
    static void write(Path testDir, Path reportDir) {
        reportDir.mkdirs()
        TestReportProcessor processor = new TestReportProcessor(new LogWrapper())
        List<TestRunResults> results = processor.scan(testDir.toFile(), 2)
        if (results.isEmpty()) {
            throw new TestReporterFailure("No XML Reports to process")
        }
        processor.write(reportDir.toFile(), new CssRagStatus(), results)
    }

    /**
     * Dummy class to ignore report generation log
     */
    static class LogWrapper implements LogFactory {

        @Override
        Log logger(String name) {
            return new NfLog()
        }

        private class NfLog implements Log {

            @Override
            void info(String msg) {}

            @Override
            void debug(String msg) {}
        }

    }
}
