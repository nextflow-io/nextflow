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

}
