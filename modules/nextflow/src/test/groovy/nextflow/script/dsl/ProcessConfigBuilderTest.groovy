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

package nextflow.script.dsl

import spock.lang.Specification
import spock.lang.Unroll

import nextflow.script.ProcessConfig
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessConfigBuilderTest extends Specification {

    def createBuilder() {
        new ProcessConfigBuilder(new ProcessConfig([:]))
    }

    @Unroll
    def 'should match selector: #SELECTOR with #TARGET' () {
        expect:
        ProcessConfigBuilder.matchesSelector(TARGET, SELECTOR) == EXPECTED

        where:
        SELECTOR        | TARGET    | EXPECTED
        'foo'           | 'foo'     | true
        'foo'           | 'bar'     | false
        '!foo'          | 'bar'     | true
        'a|b'           | 'a'       | true
        'a|b'           | 'b'       | true
        'a|b'           | 'z'       | false
        'a*'            | 'a'       | true
        'a*'            | 'aaaa'    | true
        'a*'            | 'bbbb'    | false
    }

    def 'should apply config setting for a process label' () {
        given:
        def settings = [
                'withLabel:short'  : [ cpus: 1, time: '1h'],
                'withLabel:!short' : [ cpus: 32, queue: 'cn-long'],
                'withLabel:foo'    : [ cpus: 2 ],
                'withLabel:foo|bar': [ disk: '100GB' ],
                'withLabel:gpu.+'  : [ cpus: 4 ],
        ]

        when:
        def config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithLabels(settings, ['short'])
        then:
        config.cpus == 1
        config.time == '1h'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithLabels(settings, ['long'])
        then:
        config.cpus == 32
        config.queue == 'cn-long'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithLabels(settings, ['foo'])
        then:
        config.cpus == 2
        config.disk == '100GB'
        config.queue == 'cn-long'
        config.size() == 3

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithLabels(settings, ['bar'])
        then:
        config.cpus == 32
        config.disk == '100GB'
        config.queue == 'cn-long'
        config.size() == 3

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithLabels(settings, ['gpu-1'])
        then:
        config.cpus == 4
        config.queue == 'cn-long'
        config.size() == 2

    }


    def 'should apply config setting for a process name' () {
        given:
        def settings = [
                'withName:alpha'        : [ cpus: 1, time: '1h'],
                'withName:delta'        : [ cpus: 2 ],
                'withName:delta|gamma'  : [ disk: '100GB' ],
                'withName:omega.+'      : [ cpus: 4 ],
        ]

        when:
        def config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithName(settings, 'xx')
        then:
        config.size() == 0

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithName(settings, 'alpha')
        then:
        config.cpus == 1
        config.time == '1h'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithName(settings, 'delta')
        then:
        config.cpus == 2
        config.disk == '100GB'
        config.size() == 2

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithName(settings, 'gamma')
        then:
        config.disk == '100GB'
        config.size() == 1

        when:
        config = new ProcessConfig([:])
        new ProcessConfigBuilder(config).applyConfigSelectorWithName(settings, 'omega_x')
        then:
        config.cpus == 4
        config.size() == 1
    }


    def 'should apply config process defaults' () {

        when:
        def builder = createBuilder()
        builder.queue 'cn-el6'
        builder.memory '10 GB'
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest'
        )
        def config = builder.getConfig()

        then:
        config.queue == 'cn-el6'
        config.container == 'ubuntu:latest'
        config.memory == '10 GB'
        config.cacheable == true



        when:
        builder = createBuilder()
        builder.container null
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        config = builder.getConfig()
        then:
        config.queue == 'def-queue'
        config.container == null
        config.maxRetries == 5



        when:
        builder = createBuilder()
        builder.maxRetries 10
        builder.applyConfigDefaults(
                queue: 'def-queue',
                container: 'ubuntu:latest',
                maxRetries: 5
        )
        config = builder.getConfig()
        then:
        config.queue == 'def-queue'
        config.container == 'ubuntu:latest'
        config.maxRetries == 10
    }


    def 'should apply pod configs' () {

        when:
        def builder = createBuilder()
        builder.applyConfigDefaults( pod: [secret: 'foo', mountPath: '/there'] )
        then:
        builder.getConfig().pod == [
                [secret: 'foo', mountPath: '/there']
        ]

        when:
        builder = createBuilder()
        builder.applyConfigDefaults( pod: [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ] )
        then:
        builder.getConfig().pod == [
                [secret: 'foo', mountPath: '/here'],
                [secret: 'bar', mountPath: '/there']
        ]

    }

    def 'should not apply config on negative label' () {
        given:
        def settings = [
                'withLabel:foo': [ cpus: 2 ],
                'withLabel:!foo': [ cpus: 4 ],
                'withLabel:!nodisk_.*': [ disk: '100.GB']
        ]

        when:
        def config = new ProcessConfig(label: ['foo', 'other'])
        new ProcessConfigBuilder(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 2
        config.disk == '100.GB'

        when:
        config = new ProcessConfig(label: ['foo', 'other', 'nodisk_label'])
        new ProcessConfigBuilder(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 2
        !config.disk

        when:
        config = new ProcessConfig(label: ['other', 'nodisk_label'])
        new ProcessConfigBuilder(config).applyConfig(settings, "processName", null, null)
        then:
        config.cpus == 4
        !config.disk

    }
}
