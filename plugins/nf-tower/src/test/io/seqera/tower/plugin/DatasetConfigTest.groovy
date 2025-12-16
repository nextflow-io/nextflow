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
 *
 */

package io.seqera.tower.plugin

import spock.lang.Specification

/**
 * Test DatasetConfig
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
class DatasetConfigTest extends Specification {

    def 'should create default config'() {
        when:
        def config = new DatasetConfig()

        then:
        !config.enabled
        config.createMode == 'auto'
        config.namePattern == '${workflow.runName}-outputs'
        config.perOutput.isEmpty()
    }

    def 'should create config from map'() {
        given:
        def opts = [
            enabled: true,
            createMode: 'existing',
            namePattern: 'custom-${output.name}',
            perOutput: [
                'my_output': [
                    datasetId: 'dataset-123',
                    enabled: true
                ]
            ]
        ]

        when:
        def config = new DatasetConfig(opts)

        then:
        config.enabled
        config.createMode == 'existing'
        config.namePattern == 'custom-${output.name}'
        config.perOutput.size() == 1
    }

    def 'should get output config'() {
        given:
        def opts = [
            perOutput: [
                'output1': [datasetId: 'dataset-123'],
                'output2': [enabled: false]
            ]
        ]
        def config = new DatasetConfig(opts)

        expect:
        config.getOutputConfig('output1').datasetId == 'dataset-123'
        config.getOutputConfig('output2').enabled == false
        config.getOutputConfig('output3').isEmpty()
    }

    def 'should check if enabled for output'() {
        given:
        def opts = [
            enabled: true,
            perOutput: [
                'output1': [enabled: false],
                'output2': [datasetId: 'dataset-123']
            ]
        ]
        def config = new DatasetConfig(opts)

        expect:
        !config.isEnabledForOutput('output1')  // explicitly disabled
        config.isEnabledForOutput('output2')   // enabled by default
        config.isEnabledForOutput('output3')   // enabled by default
    }

    def 'should check if disabled globally'() {
        given:
        def opts = [
            enabled: false,
            perOutput: [
                'output1': [datasetId: 'dataset-123']
            ]
        ]
        def config = new DatasetConfig(opts)

        expect:
        !config.isEnabledForOutput('output1')  // disabled globally
    }

    def 'should get dataset ID'() {
        given:
        def opts = [
            perOutput: [
                'output1': [datasetId: 'dataset-123'],
                'output2': [enabled: true]
            ]
        ]
        def config = new DatasetConfig(opts)

        expect:
        config.getDatasetId('output1') == 'dataset-123'
        config.getDatasetId('output2') == null
        config.getDatasetId('output3') == null
    }

    def 'should check auto-create mode'() {
        expect:
        new DatasetConfig([createMode: 'auto']).isAutoCreateEnabled()
        !new DatasetConfig([createMode: 'existing']).isAutoCreateEnabled()
        new DatasetConfig().isAutoCreateEnabled()  // default is 'auto'
    }

    def 'should handle empty config'() {
        when:
        def config = new DatasetConfig([:])

        then:
        !config.enabled
        config.createMode == 'auto'
        config.namePattern == '${workflow.runName}-outputs'
        config.perOutput.isEmpty()
    }

    def 'should handle null values'() {
        given:
        def opts = [
            enabled: null,
            createMode: null,
            namePattern: null,
            perOutput: null
        ]

        when:
        def config = new DatasetConfig(opts)

        then:
        !config.enabled
        config.createMode == 'auto'
        config.namePattern == '${workflow.runName}-outputs'
        config.perOutput.isEmpty()
    }

}
