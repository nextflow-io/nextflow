/*
 * Copyright 2013-2026, Seqera Labs
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

import java.nio.file.Files
import java.nio.file.Path

import ch.qos.logback.classic.Level
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.util.HistoryFile
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
class CmdLogFileTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def setup() { SysEnv.push([:]) }
    def cleanup() { SysEnv.pop() }

    private static final String ESC = ''

    private static final String SAMPLE_LOG = """\
May-26 14:30:01.000 [main] DEBUG nextflow.cli.Launcher - \$> nextflow run hello
May-26 14:30:01.100 [main] DEBUG nextflow.Session - Session UUID: aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee
May-26 14:30:01.200 [main] INFO  nextflow.Session - Starting session
May-26 14:30:02.000 [Task submitter] INFO  nextflow.Session - Task submitted
May-26 14:30:03.000 [main] WARN  nextflow.processor.TaskProcessor - Something looked odd
May-26 14:30:04.000 [main] ERROR nextflow.Session - Pipeline failed
java.lang.RuntimeException: oh no
\tat foo.Bar.baz(Bar.java:42)
\tat foo.Bar.qux(Bar.java:21)
May-26 14:30:05.000 [main] INFO  nextflow.Session - Done
""".stripIndent()

    def 'should print full file when given an explicit path'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()]).run()

        then:
        final lines = capture.toString().readLines()
        lines.any { it.contains('Session UUID') }
        lines.any { it.contains('Starting session') }
        lines.any { it.contains('Pipeline failed') }
        lines.any { it == 'java.lang.RuntimeException: oh no' }
        lines.any { it.contains('Bar.java:42') }

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should filter out DEBUG entries when level threshold is INFO'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], level: 'INFO').run()

        then:
        final out = capture.toString()
        !out.contains('DEBUG nextflow.cli.Launcher')
        !out.contains('DEBUG nextflow.Session')
        out.contains('INFO  nextflow.Session - Starting session')
        out.contains('WARN  nextflow.processor.TaskProcessor')
        out.contains('ERROR nextflow.Session - Pipeline failed')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should keep stack trace continuation lines with their parent entry'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], level: 'ERROR').run()

        then:
        final out = capture.toString()
        out.contains('ERROR nextflow.Session - Pipeline failed')
        out.contains('java.lang.RuntimeException: oh no')
        out.contains('Bar.java:42')
        out.contains('Bar.java:21')
        !out.contains('Starting session')
        !out.contains('Done')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should keep only the last N entries'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        def buf = new StringBuilder()
        (1..10).each { i ->
            buf.append(String.format("May-26 14:30:%02d.000 [main] INFO  nextflow.Session - entry %d%n", i, i))
        }
        logFile.text = buf.toString()

        when:
        new CmdLogFile(args: [logFile.toString()], lines: 3).run()

        then:
        final out = capture.toString().readLines().findAll { it.contains('entry') }
        out.size() == 3
        out[0].contains('entry 8')
        out[1].contains('entry 9')
        out[2].contains('entry 10')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should keep N entries including multi-line stack traces'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        // last 2 entries: ERROR with stack trace, then INFO Done
        new CmdLogFile(args: [logFile.toString()], lines: 2).run()

        then:
        final out = capture.toString()
        out.contains('Pipeline failed')
        out.contains('Bar.java:42')
        out.contains('Done')
        !out.contains('Starting session')
        !out.contains('Task submitted')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should resolve an older run name to a rotated .nextflow.log file by session UUID'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final uuidNew = UUID.fromString('aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee')
        final uuidOld = UUID.fromString('11111111-2222-3333-4444-555555555555')

        // current run is in .nextflow.log; previous run is rotated to .nextflow.log.1
        tmp.resolve('.nextflow.log').text = SAMPLE_LOG  // contains uuidNew
        tmp.resolve('.nextflow.log.1').text = """\
May-26 12:00:00.000 [main] DEBUG nextflow.Session - Session UUID: ${uuidOld}
May-26 12:00:01.000 [main] INFO  nextflow.Session - older run
""".stripIndent()

        final history = new HistoryFile(tmp.resolve('history').toFile())
        history.write('old_run', uuidOld, 'rev1', 'run hello')
        history.write('new_run', uuidNew, 'rev1', 'run hello')

        when:
        new CmdLogFile(args: ['old_run'], history: history, currentDir: tmp).run()

        then:
        final out = capture.toString()
        out.contains('older run')
        !out.contains('Starting session')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should resolve a current run name to the active .nextflow.log file'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final uuidNew = UUID.fromString('aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee')
        final uuidOld = UUID.fromString('11111111-2222-3333-4444-555555555555')

        tmp.resolve('.nextflow.log').text = SAMPLE_LOG
        tmp.resolve('.nextflow.log.1').text = """\
May-26 12:00:00.000 [main] DEBUG nextflow.Session - Session UUID: ${uuidOld}
May-26 12:00:01.000 [main] INFO  nextflow.Session - older run
""".stripIndent()

        final history = new HistoryFile(tmp.resolve('history').toFile())
        history.write('old_run', uuidOld, 'rev1', 'run hello')
        history.write('new_run', uuidNew, 'rev1', 'run hello')

        when:
        new CmdLogFile(args: ['new_run'], history: history, currentDir: tmp).run()

        then:
        final out = capture.toString()
        out.contains('Starting session')
        !out.contains('older run')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should error when run name is not in history'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final history = new HistoryFile(tmp.resolve('history').toFile())
        history.write('only_run', UUID.randomUUID(), 'rev1', 'run')

        when:
        new CmdLogFile(args: ['nope'], history: history, currentDir: tmp).run()

        then:
        final ex = thrown(AbortOperationException)
        ex.message.contains('Unknown run name')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should error when history is empty and no path given'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final history = new HistoryFile(tmp.resolve('history').toFile())

        when:
        new CmdLogFile(args: ['something'], history: history, currentDir: tmp).run()

        then:
        final ex = thrown(AbortOperationException)
        ex.message.contains('history is empty') || ex.message.contains('no pipeline')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should error when run name resolves but no matching log file is found'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final uuid = UUID.randomUUID()
        // history has the run, but there is no .nextflow.log file containing the UUID
        final history = new HistoryFile(tmp.resolve('history').toFile())
        history.write('orphan', uuid, 'rev1', 'run')
        tmp.resolve('.nextflow.log').text = "May-26 12:00:00.000 [main] INFO  nextflow.Session - unrelated\n"

        when:
        new CmdLogFile(args: ['orphan'], history: history, currentDir: tmp).run()

        then:
        final ex = thrown(AbortOperationException)
        ex.message.contains('Cannot find')
        ex.message.contains('orphan')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should reject too many arguments'() {
        when:
        new CmdLogFile(args: ['a', 'b']).run()

        then:
        final ex = thrown(AbortOperationException)
        ex.message.contains('Too many arguments')
    }

    def 'should reject invalid level'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], level: 'BANANA').run()

        then:
        final ex = thrown(AbortOperationException)
        ex.message.contains('Invalid log level')

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'parseLevel handles all supported levels case-insensitively'() {
        expect:
        CmdLogFile.parseLevel(input) == expected

        where:
        input    | expected
        null     | null
        ''       | null
        'trace'  | Level.TRACE
        'Debug'  | Level.DEBUG
        'INFO'   | Level.INFO
        'warn'   | Level.WARN
        'ERROR'  | Level.ERROR
    }

    // ----------------------------------------------------------------------
    // ANSI stripping
    // ----------------------------------------------------------------------

    def 'should strip ANSI codes from log content by default'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = "May-26 14:30:01.000 [main] INFO  nextflow.Session - hello ${ESC}[31mred${ESC}[0m world\n"

        when:
        new CmdLogFile(args: [logFile.toString()]).run()

        then:
        final out = capture.toString()
        out.contains('hello red world')
        !out.contains(ESC)

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'should preserve ANSI codes in content when -keep-ansi is set'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = "May-26 14:30:01.000 [main] INFO  nextflow.Session - hello ${ESC}[31mred${ESC}[0m world\n"

        when:
        new CmdLogFile(args: [logFile.toString()], keepAnsi: true).run()

        then:
        final out = capture.toString()
        out.contains("hello ${ESC}[31mred${ESC}[0m world")

        cleanup:
        tmp.toFile().deleteDir()
    }

    // ----------------------------------------------------------------------
    // ANSI output integration (per-level indicator + grammar)
    // ----------------------------------------------------------------------

    def 'should emit ANSI indicators and grammar styling when color override is true'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], colorOverride: true).run()

        then:
        final out = capture.toString()
        // Indicator blocks per level
        out.contains(LogFileFormatter.INDICATOR_INFO)             // black cell for INFO
        out.contains(LogFileFormatter.INDICATOR_WARN_FILL)        // 2 yellow cells on WARN entry-start
        out.contains(LogFileFormatter.INDICATOR_ERROR_FILL)       // 2 red cells on ERROR entry-start
        out.contains(LogFileFormatter.INDICATOR_ERROR)            // 1 red cell on ERROR continuations (stack trace)
        // TRACE/DEBUG entries get full-line dim
        out.contains("${ESC}[2m")
        // grammar styling on INFO lines: dim timestamp, cyan thread, green logger
        out.contains("${ESC}[2mMay-26 14:30:01.200${ESC}[0m")
        out.contains("${ESC}[32mnextflow.Session${ESC}[0m")
        // INFO level token stays plain
        !out.contains("${ESC}[1mINFO${ESC}[0m")
        // WARN/ERROR emit the full header region as one emphasis run
        out.contains("${ESC}[30;43mMay-26 14:30:03.000 [main] WARN${ESC}[0m")
        out.contains("${ESC}[37;41mMay-26 14:30:04.000 [main] ERROR${ESC}[0m")

        cleanup:
        tmp.toFile().deleteDir()
    }

    def '-no-ansi suppresses styling even when color override is true'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], colorOverride: true, noAnsi: true).run()

        then:
        final out = capture.toString()
        !out.contains(ESC)

        cleanup:
        tmp.toFile().deleteDir()
    }

    def 'NO_COLOR env var suppresses styling even when color override is true'() {
        given:
        SysEnv.push([NO_COLOR: '1'])
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], colorOverride: true).run()

        then:
        final out = capture.toString()
        !out.contains(ESC)

        cleanup:
        tmp.toFile().deleteDir()
        SysEnv.pop()
    }

    def 'default behaviour in test env emits no ANSI (no TTY)'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()]).run()

        then:
        final out = capture.toString()
        !out.contains(ESC)

        cleanup:
        tmp.toFile().deleteDir()
    }

    // ----------------------------------------------------------------------
    // Pager integration
    // ----------------------------------------------------------------------

    def 'resolvePagerCommand defaults to less -FR when no env var is set'() {
        expect:
        CmdLogFile.resolvePagerCommand() == ['less', '-FR']
    }

    def 'resolvePagerCommand honors PAGER env var'() {
        given:
        SysEnv.push([PAGER: 'most -s'])

        expect:
        CmdLogFile.resolvePagerCommand() == ['most', '-s']
    }

    def 'resolvePagerCommand prefers NXF_PAGER over PAGER'() {
        given:
        SysEnv.push([NXF_PAGER: 'bat --paging always', PAGER: 'less'])

        expect:
        CmdLogFile.resolvePagerCommand() == ['bat', '--paging', 'always']
    }

    def 'parsePagerCommand appends -FR to a bare less invocation'() {
        expect:
        CmdLogFile.parsePagerCommand('less')          == ['less', '-FR']
        CmdLogFile.parsePagerCommand('/usr/bin/less') == ['/usr/bin/less', '-FR']
    }

    def 'parsePagerCommand passes through less with explicit args'() {
        expect:
        CmdLogFile.parsePagerCommand('less -R')   == ['less', '-R']
        CmdLogFile.parsePagerCommand('less -RFX') == ['less', '-RFX']
    }

    def 'parsePagerCommand handles arbitrary pagers without augmentation'() {
        expect:
        CmdLogFile.parsePagerCommand('most')        == ['most']
        CmdLogFile.parsePagerCommand('cat')         == ['cat']
        CmdLogFile.parsePagerCommand('bat --plain') == ['bat', '--plain']
    }

    def 'parsePagerCommand returns null for empty input'() {
        expect:
        CmdLogFile.parsePagerCommand('')      == null
        CmdLogFile.parsePagerCommand('   ')   == null
    }

    def 'shouldUsePager is false when -no-pager is set'() {
        expect:
        !new CmdLogFile(noPager: true).shouldUsePager()
    }

    def 'shouldUsePager is false in follow mode'() {
        expect:
        !new CmdLogFile(follow: true).shouldUsePager()
    }

    def 'shouldUsePager is false when stdout is not a terminal (no TTY in test env)'() {
        expect:
        !new CmdLogFile().shouldUsePager()
    }

    def 'should preserve indicator + grammar on the last-N path'() {
        given:
        final tmp = Files.createTempDirectory('cmdlogfile-')
        final logFile = tmp.resolve('.nextflow.log')
        logFile.text = SAMPLE_LOG

        when:
        new CmdLogFile(args: [logFile.toString()], lines: 2, colorOverride: true).run()

        then:
        final out = capture.toString()
        // last two entries: ERROR (with stack trace) + INFO Done
        out.contains(LogFileFormatter.INDICATOR_ERROR_FILL)   // ERROR entry-start
        out.contains(LogFileFormatter.INDICATOR_ERROR)        // ERROR continuation (stack trace)
        out.contains(LogFileFormatter.INDICATOR_INFO + ESC + '[2mMay-26 14:30:05.000')
        // exception line highlighted under the ERROR entry
        out.contains("${ESC}[1;31mjava.lang.RuntimeException${ESC}[0m")

        cleanup:
        tmp.toFile().deleteDir()
    }
}
