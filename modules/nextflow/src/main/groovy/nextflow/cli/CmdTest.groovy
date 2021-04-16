/*
 * Copyright 2020-2021, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.cli
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.NextflowMeta
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.script.ScriptRunner
import nextflow.script.testflow.HtmlRenderer
import nextflow.script.testflow.TestRun
import nextflow.script.testflow.TestSuite
import org.fusesource.jansi.Ansi

import java.nio.file.FileSystems
import java.nio.file.Path
import java.nio.file.Paths
import java.time.Duration
import java.time.Instant

/**
 * CLI sub-command TEST
 *
 * @author Jordi Deu-Pons <jordi@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Run all the project tests")
class CmdTest extends CmdBase implements HubOptions {

    static final public NAME = 'test'

    @Parameter(description = 'Project name or repository url')
    List<String> args

    @Parameter(names=['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names=['-offline'], description = 'Do not check for remote project updates')
    boolean offline = System.getenv('NXF_OFFLINE') as boolean

    @Parameter(names=['-w', '-work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @Parameter(names='-pattern', description = 'Pattern to find test scripts')
    String pattern = "**/*_test.nf"

    @Parameter(names=['-skip-report'], description = 'Do not create HTML test reports')
    boolean skipReport

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        log.info "N E X T F L O W  ~  version ${Const.APP_VER}"

        checkValidParams()

        // Test are only supported with DSL2
        NextflowMeta.instance.enableDsl2()

        final repoName = args ? args[0] : Paths.get("").toAbsolutePath().toString()
        final localPath = getLocalPath(repoName)

        List<Path> testScripts = new ArrayList<>()
        final matcher = FileSystems.getDefault().getPathMatcher("glob:${pattern}")
        localPath.eachFileRecurse(groovy.io.FileType.FILES) {
            if (matcher.matches(it)) {
                testScripts.add(it)
            }
        }

        // create the config object
        final builder = new ConfigBuilder()
                .setOptions(launcher.options)
                .setBaseDir(localPath)
        final config = builder.build()

        // -- sets the working directory
        if( workDir )
            config.workDir = workDir
        else if( !config.workDir )
            config.workDir = System.getenv().get('NXF_WORK') ?: 'work'

        final configWorkDir = config.workDir as Path
        final testResults = configWorkDir.resolve("test-results")
        TestRun testRun = new TestRun()
        for (script in testScripts) {
            String testName = localPath.relativize(script).toString()
                    .replace(FileSystems.getDefault().separator, ".")
                    .replace(".nf", "")
            log.info ""
            log.info printAnsi(Ansi.Color.WHITE,"Running test > ${testName}")
            config.workDir = configWorkDir.resolve("test-runs/${testName}")
            testRun.addSuite(
                runTest(testName, script, testResults, config)
            )
        }

        String failureMessage
        if (!testRun.empty) {

            if (testRun.failed > 0) {
                for (suite in testRun.suites) {
                    for (testcase in suite.testcase) {
                        if (testcase.failure) {
                            log.info ""
                            log.info printAnsi(Ansi.Color.RED, "Failure at ${suite.name} > ${testcase.name}")
                            log.info testcase.failure.content
                            log.info ""
                        }
                    }
                }

                failureMessage = printAnsi(Ansi.Color.RED, "TESTS FAILED")
            }

            log.info ""
            log.info printAnsi(Ansi.Color.YELLOW, "${testRun.completed} tests completed, ${testRun.skipped} skipped, ${testRun.failed} failed")
            log.info ""
        }

        // Create HTML reports
        if (!skipReport) {
            final reportDir = configWorkDir.resolve("test-reports")
            HtmlRenderer.write(testRun, reportDir)

            final reportIndex = reportDir.toAbsolutePath().resolve("index.html")
            log.info "Open test reports at: " + printAnsi(Ansi.Color.WHITE, "file://${reportIndex}")
            log.info ""
        }

        if (failureMessage) {
            throw new AbortOperationException(failureMessage)
        } else {
            log.info printAnsi(Ansi.Color.GREEN, "TESTS DONE")
        }

    }

    protected void checkValidParams() {
        if (offline && latest)
            throw new AbortOperationException("Command line options `-latest` and `-offline` cannot be specified at the same time")
    }

    private String printAnsi(Ansi.Color color, String message) {
        if (launcher.options.ansiLog)
            return Ansi.ansi().fg(color).a(message).reset()
        return message
    }

    protected TestSuite runTest(String testName, Path testScript, Path testResults, Map config) {

        Plugins.setup( config )

        final runStart = Instant.now()
        final runner = new ScriptRunner(config)
        runner.setScript(testScript)
        runner.session.profile = profile
        runner.session.commandLine = launcher.cliString
        runner.session.ansiLog = launcher.options.ansiLog
        runner.session.testName = testName
        runner.session.testResults = testResults

        runner.testFlow()

        final suite = runner.getTestSuite()
        suite.time = Duration.between(runStart, Instant.now())
        return suite
    }

    protected Path getLocalPath(String repoName) {

        /*
         * try to find a local pipeline
         */
        def baseDir = new File(repoName)
        if ( baseDir.exists() ) {
            if ( baseDir.isFile() ) {
                return baseDir.getParentFile().toPath()
            }
            if ( baseDir.isDirectory() ) {
                return baseDir.toPath()
            }
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new AssetManager(repoName, this)
        def repo = manager.getProject()

        boolean checkForUpdate = true
        if( !manager.isRunnable() || latest ) {
            if( offline )
                throw new AbortOperationException("Unknown project `$repo` -- NOTE: automatic download from remote repositories is disabled")
            log.info "Pulling $repo ..."
            def result = manager.download()
            if( result )
                log.info " $result"
            checkForUpdate = false
        }
        // checkout requested revision
        try {
            manager.checkout(revision)
            manager.updateModules()
            if( checkForUpdate && !offline )
                manager.checkRemoteStatus(manager.getCurrentRevisionAndName())
            // return the local folder
            return manager.getLocalPath().toPath()
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Unknown error accessing project `$repo` -- Repository may be corrupted: ${manager.localPath}", e)
        }
    }
}
