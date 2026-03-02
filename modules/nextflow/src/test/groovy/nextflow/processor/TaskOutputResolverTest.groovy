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

package nextflow.processor

import java.nio.file.Path

import nextflow.exception.IllegalArityException
import nextflow.exception.MissingFileException
import nextflow.exception.MissingValueException
import nextflow.script.params.v2.ProcessFileOutput
import spock.lang.Specification
import spock.lang.TempDir

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TaskOutputResolverTest extends Specification {

    @TempDir
    Path tempDir

    def makeTask(Map holder = [:]) {
        return Spy(new TaskRun(
            name: 'test_task',
            workDir: tempDir,
            config: new TaskConfig(),
            context: new TaskContext(holder: holder)
        ))
    }

    def makeFileOutput(String pattern) {
        return Mock(ProcessFileOutput) {
            getFilePattern(_) >> pattern
        }
    }

    def 'should get environment variable'() {
        given:
        def task = makeTask()
        task.getOutputEvals() >> [:]
        and:
        def resolver = new TaskOutputResolver([:], task)
        and:
        def envFile = tempDir.resolve(TaskRun.CMD_ENV)
        envFile.text = '''
            FOO=bar
            /FOO/
            BAZ=qux
            /BAZ/
            '''.stripIndent()

        when:
        def result = resolver._env('FOO')
        then:
        result == 'bar'

        when:
        resolver._env('MISSING')
        then:
        thrown(MissingValueException)
    }

    def 'should get eval result'() {
        given:
        def task = makeTask()
        task.getOutputEvals() >> [RESULT: 'echo hello']
        and:
        def resolver = new TaskOutputResolver([:], task)
        and:
        def envFile = tempDir.resolve(TaskRun.CMD_ENV)
        envFile.text = 'RESULT=hello\n'

        when:
        def result = resolver.eval('RESULT')
        then:
        result == 'hello'

        when:
        resolver.eval('MISSING')
        then:
        thrown(MissingValueException)
    }

    def 'should get single file'() {
        given:
        def declaredFiles = [out: makeFileOutput('output.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)
        and:
        def outputFile = tempDir.resolve('output.txt')
        outputFile.text = ''

        when:
        def result = resolver._file('out')
        then:
        result == outputFile
        task.outputFiles == [outputFile].toSet()
    }

    def 'should fail when file pattern is absolute path'() {
        given:
        def declaredFiles = [out: makeFileOutput('/absolute/path/output.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        resolver._file('out')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('is an absolute path')
    }

    def 'should fail when single file yields no matches'() {
        given:
        def declaredFiles = [out: makeFileOutput('missing.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        resolver._file('out')
        then:
        def e = thrown(MissingFileException)
        e.message.contains('Missing output file(s) `missing.txt`')
    }

    def 'should return null when optional single file yields no matches'() {
        given:
        def declaredFiles = [out: makeFileOutput('missing.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        def result = resolver._file([optional: true], 'out')
        then:
        result == null
    }

    def 'should fail when single file yields multiple matches'() {
        given:
        def declaredFiles = [out: makeFileOutput('*.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)
        and:
        tempDir.resolve('file1.txt').text = ''
        tempDir.resolve('file2.txt').text = ''

        when:
        resolver._file('out')
        then:
        def e = thrown(IllegalArityException)
        e.message.contains('yielded 2 files but expected only one')
    }

    def 'should get multiple files'() {
        given:
        def declaredFiles = [out: makeFileOutput('*.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)
        and:
        def file1 = tempDir.resolve('file1.txt')
        def file2 = tempDir.resolve('file2.txt')
        file1.text = ''
        file2.text = ''

        when:
        def result = resolver._files('out')
        then:
        result == [file1, file2] as Set
        task.outputFiles.containsAll(result)
    }

    def 'should fail when files pattern is absolute path'() {
        given:
        def declaredFiles = [out: makeFileOutput('/absolute/path/*.txt')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        resolver._files('out')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('is an absolute path')
    }

    def 'should fail when files glob yields no matches'() {
        given:
        def declaredFiles = [out: makeFileOutput('*.missing')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        resolver._files('out')
        then:
        def e = thrown(MissingFileException)
        e.message.contains('Missing output file(s) `*.missing`')
    }

    def 'should return empty set when optional files glob yields no matches'() {
        given:
        def declaredFiles = [out: makeFileOutput('*.missing')]
        def task = makeTask()
        def resolver = new TaskOutputResolver(declaredFiles, task)

        when:
        def result = resolver._files([optional: true], 'out')
        then:
        result == [] as Set
    }

    def 'should get stdout from file'() {
        given:
        def stdoutFile = tempDir.resolve('.command.out')
        stdoutFile.text = 'hello world'
        and:
        def task = makeTask()
        task.stdout = stdoutFile
        def resolver = new TaskOutputResolver([:], task)

        expect:
        resolver.stdout() == 'hello world'
    }

    def 'should get stdout from text'() {
        given:
        def task = makeTask()
        task.stdout = 'direct value'
        def resolver = new TaskOutputResolver([:], task)

        expect:
        resolver.stdout() == 'direct value'
    }

    def 'should fail when stdout is missing'() {
        given:
        def task = makeTask()
        task.stdout = null
        def resolver = new TaskOutputResolver([:], task)

        when:
        resolver.stdout()
        then:
        thrown(IllegalArgumentException)
    }

    def 'should fail when stdout file does not exist'() {
        given:
        def missingFile = tempDir.resolve('missing.out')
        def task = makeTask()
        task.stdout = missingFile
        def resolver = new TaskOutputResolver([:], task)

        when:
        resolver.stdout()
        then:
        thrown(Exception)
    }

    def 'should get variable from context'() {
        given:
        def task = makeTask([foo: 'bar', baz: 123])
        def resolver = new TaskOutputResolver([:], task)

        expect:
        resolver.get('foo') == 'bar'
        resolver.get('baz') == 123
    }

    def 'should fail when variable is missing from context'() {
        given:
        def task = makeTask([foo: 'bar'])
        def resolver = new TaskOutputResolver([:], task)

        when:
        resolver.get('missing')
        then:
        def e = thrown(MissingValueException)
        e.message.contains('Missing variable in process output')
    }
}
