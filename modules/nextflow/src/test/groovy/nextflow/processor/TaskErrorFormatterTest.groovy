/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.processor

import java.nio.file.Path

import nextflow.exception.FailedGuardException
import nextflow.exception.ProcessEvalException
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskErrorFormatterTest extends Specification {

    def 'should format command error'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
        }
        def error = new ProcessEvalException(
            'Command failed',
            '  echo "hello"\n  exit 1\n',
            'hello\nerror output',
            1
        )

        when:
        def result = formatter.formatCommandError([], error, task).join('\n')

        then:
        result.contains('Caused by:')
        result.contains('Command executed:')
        result.contains('echo "hello"')
        result.contains('exit 1')
        result.contains('Command exit status:')
        result.contains('1')
        result.contains('Command output:')
        result.contains('hello')
        result.contains('error output')
        result.contains('Work dir:')
        result.contains('/work/dir')
    }

    def 'should format command error with empty output'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
        }
        def error = new ProcessEvalException(
            'Command failed',
            'exit 1',
            '',
            1
        )

        when:
        def result = formatter.formatCommandError([], error, task).join('\n')

        then:
        result.contains('Command output:')
        result.contains('(empty)')
    }

    def 'should format command error without work dir'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> null
        }
        def error = new ProcessEvalException(
            'Command failed',
            'exit 1',
            'output',
            1
        )

        when:
        def result = formatter.formatCommandError([], error, task).join('\n')

        then:
        !result.contains('Work dir:')
    }

    def 'should format guard error'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
        }
        def error = new FailedGuardException(
            'Guard failed',
            'when: false',
            new IllegalStateException('Invalid condition')
        )

        when:
        def result = formatter.formatGuardError([], error, task).join('\n')

        then:
        result.contains('Caused by:')
        result.contains('When block:')
        result.contains('when:')
        result.contains('false')
        result.contains('Work dir:')
        result.contains('/work/dir')
    }

    def 'should format guard error without work dir'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> null
        }
        def error = new FailedGuardException(
            'Guard failed',
            'when: false',
            new IllegalStateException('Invalid condition')
        )

        when:
        def result = formatter.formatGuardError([], error, task).join('\n')

        then:
        !result.contains('Work dir:')
    }

    def 'should format task error with script block'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> '  echo "test"\n  exit 1\n'
            getTemplate() >> null
            getExitStatus() >> 1
            dumpStdout(_) >> ['stdout line 1', 'stdout line 2']
            dumpStderr(_) >> ['stderr line 1']
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Caused by:')
        result.contains('Command executed:')
        result.contains('echo "test"')
        result.contains('exit 1')
        result.contains('Command exit status:')
        result.contains('1')
        result.contains('Command output:')
        result.contains('stdout line 1')
        result.contains('stdout line 2')
        result.contains('Command error:')
        result.contains('stderr line 1')
        result.contains('Work dir:')
        result.contains('/work/dir')
    }

    def 'should format task error with template'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> 'echo "test"'
            getTemplate() >> Path.of('template.sh')
            getExitStatus() >> 1
            dumpStdout(_) >> ['output']
            dumpStderr(_) >> []
            dumpLogFile(_) >> ['log line 1']
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Command executed [template.sh]:')
        result.contains('Command log:')
        result.contains('log line 1')
    }

    def 'should format task error with empty stdout'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> 'exit 1'
            getTemplate() >> null
            getExitStatus() >> 1
            dumpStdout(_) >> []
            dumpStderr(_) >> ['error']
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Command output:')
        result.contains('(empty)')
    }

    def 'should format task error with exec block'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> null
            getSource() >> '  def x = 1\n  println x'
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Caused by:')
        result.contains('Source block:')
        result.contains('def x = 1')
        result.contains('println x')
        result.contains('Work dir:')
    }

    def 'should format task error with container'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> 'echo test'
            getTemplate() >> null
            getExitStatus() >> 0
            dumpStdout(_) >> ['output']
            dumpStderr(_) >> []
            isContainerEnabled() >> true
            getContainer() >> 'ubuntu:20.04'
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Container:')
        result.contains('ubuntu:20.04')
    }

    def 'should format task error without work dir'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> null
            getScript() >> 'echo test'
            getTemplate() >> null
            getExitStatus() >> 0
            dumpStdout(_) >> ['output']
            dumpStderr(_) >> []
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        !result.contains('Work dir:')
    }

    def 'should format task error with MAX_VALUE exit status'() {
        given:
        def formatter = new TaskErrorFormatter()
        def task = Mock(TaskRun) {
            getWorkDir() >> Path.of('/work/dir')
            getWorkDirStr() >> '/work/dir'
            getScript() >> 'echo test'
            getTemplate() >> null
            getExitStatus() >> Integer.MAX_VALUE
            dumpStdout(_) >> []
            dumpStderr(_) >> []
            isContainerEnabled() >> false
        }
        def error = new RuntimeException('Task failed')

        when:
        def result = formatter.formatTaskError([], error, task).join('\n')

        then:
        result.contains('Command exit status:')
        result.contains('-')
    }

}
