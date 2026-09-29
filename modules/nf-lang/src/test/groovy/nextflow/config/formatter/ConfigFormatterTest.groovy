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

package nextflow.config.formatter

import nextflow.config.control.ConfigParser
import nextflow.script.formatter.FormattingOptions
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigFormatterTest extends Specification {

    @Shared
    ConfigParser parser

    def setupSpec() {
        parser = new ConfigParser()
    }

    String format(String contents) {
        def source = parser.parse('main.nf', contents)
        assert !source.getErrorCollector().hasErrors()
        def formatter = new ConfigFormattingVisitor(source, new FormattingOptions(4, true))
        formatter.visit()
        assert formatter.getMissingComments().isEmpty()
        return formatter.toString()
    }

    boolean checkFormat(String source) {
        source = source.stripIndent()
        assert format(source) == source
        return true
    }

    boolean checkFormat(String input, String output) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert format(input) == output
        assert format(output) == output
        return true
    }

    String formatWrapped(String contents, int lineLength) {
        def source = parser.parse('main.nf', contents)
        assert !source.getErrorCollector().hasErrors()
        def formatter = new ConfigFormattingVisitor(source, new FormattingOptions(4, true, false, false, false, lineLength))
        formatter.visit()
        assert formatter.getMissingComments().isEmpty()
        return formatter.toString()
    }

    boolean checkFormatWrapped(String input, String output, int lineLength = 40) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert formatWrapped(input, lineLength) == output
        assert formatWrapped(output, lineLength) == output
        return true
    }

    def 'should format a config assignment' () {
        expect:
        checkFormat(
            '''\
            process . clusterOptions={"--cpus ${task.cpus}"}
            ''',
            '''\
            process.clusterOptions = { "--cpus ${task.cpus}" }
            '''
        )
    }

    def 'should format a config block' () {
        expect:
        checkFormat(
            '''\
            process{clusterOptions={"--cpus ${task.cpus}"}}
            ''',
            '''\
            process {
                clusterOptions = { "--cpus ${task.cpus}" }
            }
            '''
        )
    }

    def 'should format a config apply block' () {
        expect:
        checkFormat(
            '''\
            plugins{id'nf-hello'}
            ''',
            '''\
            plugins {
                id 'nf-hello'
            }
            '''
        )
    }

    /// COMMENTS

    def 'should preserve comments in a config file' () {
        expect:
        checkFormat(
            '''\
            // about the process scope
            process {
                cpus = 2 // how many
                // dangling
            }

            // about the plugins
            plugins {
                id 'nf-hello' // a plugin
                // dangling
            }

            // trailing comment at the end of the file
            '''
        )
    }

    def 'should preserve a comment on an include' () {
        expect:
        checkFormat(
            '''\
            // about the include
            includeConfig 'other.config' // over here
            '''
        )
    }

    def 'should re-indent a multi-line string with its statement' () {
        expect:
        checkFormat(
            '''\
            process {
              beforeScript = """
              echo one
              """
            }
            ''',
            '''\
            process {
                beforeScript = """
                echo one
                """
            }
            '''
        )
    }

    // -- line-length wrapping (issue #26)

    def 'should wrap a long config assignment and leave a short one alone' () {
        expect:
        checkFormatWrapped(
            '''\
            process.ext.args = [alpha, beta, gamma, delta, epsilon, zeta]
            ''',
            '''\
            process.ext.args = [
                alpha,
                beta,
                gamma,
                delta,
                epsilon,
                zeta,
            ]
            '''
        )
        checkFormatWrapped(
            '''\
            process.ext.args = [alpha, beta]
            ''',
            '''\
            process.ext.args = [alpha, beta]
            '''
        )
    }

    def 'should not wrap a long config assignment when maxLineLength is disabled' () {
        expect:
        checkFormatWrapped(
            '''\
            process.ext.args = [alpha, beta, gamma, delta, epsilon, zeta]
            ''',
            '''\
            process.ext.args = [alpha, beta, gamma, delta, epsilon, zeta]
            ''',
            0
        )
    }

    // -- fmt: skip / fmt: off / fmt: on directives (issue #75)

    def 'should keep a config assignment verbatim with fmt: skip' () {
        expect:
        checkFormat(
            '''\
            process.cpus=2  // fmt: skip
            process.memory='4 GB'
            ''',
            '''\
            process.cpus=2  // fmt: skip
            process.memory = '4 GB'
            '''
        )
    }

    def 'should keep a region of a config file verbatim between fmt: off and fmt: on' () {
        expect:
        checkFormat(
            '''\
            process.cpus = 2

            // fmt: off
            docker{
            enabled=true
            }
            // fmt: on

            process.memory = '4 GB'
            '''
        )
    }

}
