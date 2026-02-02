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

    def 'should warn about invalid config options' () {
        when:
        new ConfigValidator().validate([
            wokDir: 'work',
            workDir: 'work',
            process: [
                cpu: 2,
                cpus: 2
            ]
        ])
        then:
        capture.toString().contains('Unrecognized config option \'wokDir\'')
        capture.toString().contains('Unrecognized config option \'process.cpu\'')
        !capture.toString().contains('Unrecognized config option \'workDir\'')
        !capture.toString().contains('Unrecognized config option \'process.cpus\'')
    }

    def 'should validate config options in profiles' () {
        when:
        new ConfigValidator().validate([
            profiles: [
                test: [
                    workDir: 'work',
                    process: [
                        cpus: 2
                    ]
                ]
            ]
        ])
        then:
        !capture.toString().contains('Unrecognized config option')

        when:
        new ConfigValidator().validate([
            profiles: [
                test: [
                    wokDir2: 'work',
                    process: [
                        cpu2: 2
                    ]
                ]
            ]
        ])
        then:
        capture.toString().contains('Unrecognized config option \'wokDir2\'')
        capture.toString().contains('Unrecognized config option \'process.cpu2\'')
    }

    def 'should warn about invalid env config options' () {
        when:
        new ConfigValidator().validate([
            env: [
                FOO: '/something',
                NXF_ANSI_SUMMARY: 'true',
                NXF_DEBUG: 'true'
            ]
        ])

        then:
        capture.toString().contains('the following environment variable in the config will be ignored: \'NXF_ANSI_SUMMARY\'')
        !capture.toString().contains('the following environment variable in the config will be ignored: \'FOO\'')
        !capture.toString().contains('the following environment variable in the config will be ignored: \'NXF_DEBUG\'')
    }

    def 'should ignore process selectors' () {
        when:
        new ConfigValidator().validate([
            process: [
                'withLabel:foobar': [
                    cpus: 2
                ],
                'withName:foobar': [
                    cpus: 2
                ],
                "withName:'.*TASK.*'": [
                    cpus: 2
                ]
            ]
        ])
        then:
        !capture.toString().contains('Unrecognized config option')
    }

    def 'should support map options' () {
        when:
        new ConfigValidator().validate([
            k8s: [
                pod: [env: 'MESSAGE', value: 'hello world']
            ]
        ])
        then:
        !capture.toString().contains('Unrecognized config option')

        when:
        new ConfigValidator().validate([
            process: [
                publishDir: [
                    path: { "results/foo" },
                    mode: 'copy',
                    saveAs: { filename -> filename }
                ]
            ]
        ])
        then:
        !capture.toString().contains('Unrecognized config option')

        when:
        new ConfigValidator().validate([
            process: [
                resourceLimits: [
                    cpus: 4,
                    memory: '10GB',
                    time: '1.h'
                ]
            ]
        ])
        then:
        !capture.toString().contains('Unrecognized config option')
    }

    def 'should not validate core plugin config when plugin is not loaded' () {
        when:
        new ConfigValidator().validate([
            cloudcache: [
                enabled: true,
                path: 's3://bucket/cache'
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'cloudcache'")
    }

}
