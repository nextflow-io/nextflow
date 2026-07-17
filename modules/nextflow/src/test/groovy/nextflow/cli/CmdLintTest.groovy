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

import nextflow.exception.AbortOperationException
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.SourceUnit
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CmdLintTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    CmdLint createCmdLint(Map opts = null) {
        def cmd = new CmdLint()
        cmd.launcher = new Launcher(
            options: new CliOptions(opts ?: [ansiLog: false])
        )
        return cmd
    }

    SourceUnit fileSource(Path path) {
        return new SourceUnit(path.toFile(), new CompilerConfiguration(), null, null)
    }

    def 'should read files from positional args and -files-from option' () {

        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('main.nf').text = ''
        dir.resolve('nextflow.config').text = ''
        and:
        def fileList = dir.resolve('files.txt')
        fileList.text = dir.resolve('nextflow.config').toString() + '\n'

        when:
        def cmd = createCmdLint()
        cmd.args = [dir.resolve('main.nf').toString()]
        cmd.filesFrom = fileList.toFile().toString()
        cmd.run()

        then:
        capture.toString().contains(dir.resolve('main.nf').toString())
        capture.toString().contains(dir.resolve('nextflow.config').toString())

        cleanup:
        dir?.deleteDir()
    }

    def 'should throw error when -files-from file is empty' () {

        given:
        def dir = Files.createTempDirectory('test')
        def fileList = dir.resolve('files.txt')
        fileList.text = '\n'

        when:
        def cmd = createCmdLint()
        cmd.filesFrom = fileList.toFile().toString()
        cmd.run()

        then:
        thrown(AbortOperationException)

        cleanup:
        dir?.deleteDir()
    }

    def 'should report compilation errors' () {

        given:
        def dir = Files.createTempDirectory('test')

        dir.resolve('main.nf').text = '''\
            process HELLO {

                script:
                """
                ${
                    params.is_paired_end
                        ? "..."
                        : "..."
                }
                """
            }
            '''

        dir.resolve('nextflow.config').text = '''\
            process {
                withLabel:
                'bambino' {
                    container = "..."
                }
            }
            '''

        when:
        def cmd = createCmdLint()
        cmd.args = ['.']
        cmd.root = dir
        cmd.run()

        then:
        thrown(AbortOperationException)
        and:
        capture.toString().contains("Error main.nf:5:19: Unexpected input: '\\n'")
        capture.toString().contains("│   5 |                 \${")
        capture.toString().contains("╰     |                   ^")
        capture.toString().contains("Error nextflow.config:2:27: Unexpected input: '\\n'")
        capture.toString().contains("│   2 |                 withLabel:")
        capture.toString().contains("╰     |                           ^")

        cleanup:
        dir?.deleteDir()
    }

    def 'should add lib directory to class loader' () {

        given:
        def dir = Files.createTempDirectory('test')

        dir.resolve('main.nf').text = '''\
            println Utils.hello()
            '''

        dir.resolve('lib').mkdir()
        dir.resolve('lib/Utils.groovy').text = '''\
            class Utils {

                String hello() {
                    return 'Hello!'
                }
            }
            '''

        when:
        def cmd = createCmdLint()
        cmd.args = [dir.toString()]
        cmd.projectDir = dir.toString()
        cmd.run()

        then:
        noExceptionThrown()

        cleanup:
        dir?.deleteDir()
    }

    def 'should suppress progress with -q flag'() {

        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('main.nf').text = ''

        when:
        def cmd = createCmdLint(ansiLog: false, quiet: true)
        cmd.args = [dir.toFile().toString()]
        cmd.run()

        then:
        noExceptionThrown()
        and:
        !capture.toString().contains("Linting Nextflow code")
        !capture.toString().contains("Linting:")
        and:
        capture.toString().contains("Nextflow linting complete")

        cleanup:
        dir?.deleteDir()
    }

    def 'should still show errors when progress is suppressed' () {

        given:
        def dir = Files.createTempDirectory('test')

        dir.resolve('main.nf').text = '''\
            printx 'hello'
            '''

        when:
        def cmd = createCmdLint(ansiLog: false, quiet: true)
        cmd.args = [dir.toFile().toString()]
        cmd.run()

        then:
        thrown(AbortOperationException)
        and:
        !capture.toString().contains("Linting Nextflow code")
        !capture.toString().contains("Linting:")
        and:
        capture.toString().contains("Error")
        capture.toString().contains("`printx` is not defined")
        capture.toString().contains("Nextflow linting complete")

        cleanup:
        dir?.deleteDir()
    }

    def 'should exclude sources under an excluded directory' () {

        given:
        def dir = Files.createTempDirectory('test')
        def included = dir.resolve('modules/local/foo/main.nf')
        def excluded = dir.resolve('modules/nf-core/bar/main.nf')
        Files.createDirectories(included.parent)
        Files.createDirectories(excluded.parent)
        included.text = ''
        excluded.text = ''
        and:
        def cmd = createCmdLint()
        cmd.excludePatterns = ['nf-core']

        expect:
        // a file directly under a non-excluded directory is kept
        !cmd.isExcludedSource(fileSource(included))
        // a file reached indirectly (e.g. via include) under an excluded
        // ancestor directory is excluded, even with an absolute path
        cmd.isExcludedSource(fileSource(excluded))

        cleanup:
        dir?.deleteDir()
    }

    def 'should report indirectly-included sources with a relative path' () {

        given:
        def cwd = Path.of('.').toAbsolutePath()
        // create the file under the working directory so it can be relativized
        def dir = Files.createTempDirectory(cwd, 'test')
        def file = dir.resolve('modules/nf-core/bar/main.nf')
        Files.createDirectories(file.parent)
        file.text = ''
        and:
        def cmd = createCmdLint()
        def source = fileSource(file)

        expect:
        // the source is backed by an absolute file, as if resolved from an include
        source.getName() == file.toString()
        and:
        // but it is reported relative to the working directory
        cmd.relativeName(source) == cwd.relativize(file).toString()

        cleanup:
        dir?.deleteDir()
    }

}
