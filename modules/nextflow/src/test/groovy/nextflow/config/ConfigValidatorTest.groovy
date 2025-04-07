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

package nextflow.config

import nextflow.SysEnv
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ConfigValidatorTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def setupSpec() {
        SysEnv.push(NXF_SYNTAX_PARSER: 'v2')
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should warn about invalid config options' () {
        given:
        def config = new ConfigMap([
            wokDir: 'work',
            workDir: 'work',
            process: [
                cpu: 2,
                cpus: 2
            ]
        ])

        when:
        new ConfigValidator().validate(config)
        then:
        capture.toString().contains('Unrecognized config option \'wokDir\'')
        capture.toString().contains('Unrecognized config option \'process.cpu\'')
        !capture.toString().contains('Unrecognized config option \'workDir\'')
        !capture.toString().contains('Unrecognized config option \'process.cpus\'')
    }

    def 'should warn about invalid env config options' () {
        when:
        new ConfigValidator().validate(new ConfigMap([
            env: [
                FOO: '/something',
                NXF_ANSI_SUMMARY: 'true',
                NXF_DEBUG: 'true'
            ]
        ]))

        then:
        capture.toString().contains('the following environment variable in the config will be ignored: \'NXF_ANSI_SUMMARY\'')
        !capture.toString().contains('the following environment variable in the config will be ignored: \'FOO\'')
        !capture.toString().contains('the following environment variable in the config will be ignored: \'NXF_DEBUG\'')
    }

}
