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

    String format(FormattingOptions options, String contents) {
        def source = parser.parse('main.nf', contents)
        assert !source.getErrorCollector().hasErrors()
        def formatter = new ConfigFormattingVisitor(source, options, contents)
        formatter.visit()
        return formatter.toString()
    }

    String format(String contents) {
        return format(new FormattingOptions(4, true), contents)
    }

    boolean checkFormat(FormattingOptions options, String input, String output) {
        input = input.stripIndent()
        output = output.stripIndent()
        assert format(options, input) == output
        assert format(options, output) == output
        return true
    }

    boolean checkFormat(String input, String output) {
        return checkFormat(new FormattingOptions(4, true), input, output)
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

    def 'should preserve comments when formatting a config file' () {
        expect:
        checkFormat(
            '''\
            // header comment
            process {
                cpus = 4 // inline comment
                // dangling in block
            }

            // trailing comment at EOF
            // process { memory = '8GB' }
            ''',
            '''\
            // header comment
            process {
                cpus = 4 // inline comment
                // dangling in block
            }

            // trailing comment at EOF
            // process { memory = '8GB' }
            '''
        )
        checkFormat(
            '''\
            profiles {
                // the test profile
                test {
                    params.x = 1
                    // done
                }
                // dangling in profiles
            }

            process {
                // nothing here yet
            }
            ''',
            '''\
            profiles {
                // the test profile
                test {
                    params.x = 1
                    // done
                }
                // dangling in profiles
            }

            process {
                // nothing here yet
            }
            '''
        )
        // CRLF line endings
        checkFormat(
            "process {\r\n    // why\r\n    cpus = 4\r\n}\r\n",
            '''\
            process {
                // why
                cpus = 4
            }
            '''
        )
    }

}
