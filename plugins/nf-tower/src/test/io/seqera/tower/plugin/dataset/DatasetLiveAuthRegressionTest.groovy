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

package io.seqera.tower.plugin.dataset

import nextflow.Global
import nextflow.Session
import spock.lang.PendingFeature
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Live regression for dataset file reads backed by Seqera API URLs.
 *
 * Requires env vars:
 * - TOWER_ACCESS_TOKEN
 * - TOWER_WORKSPACE_ID
 * - DATASET_LIVE_NAMES (comma-separated)
 */
@Requires({ env['TOWER_ACCESS_TOKEN'] && env['TOWER_WORKSPACE_ID'] && env['DATASET_LIVE_NAMES'] })
class DatasetLiveAuthRegressionTest extends Specification {

    def cleanup() {
        Global.session = null
    }

    @PendingFeature(reason = 'Dataset provider does not yet forward Tower bearer auth to resolved HTTP dataset URLs')
    @Unroll
    def 'should read live dataset via provider using bearer auth - #datasetName'() {
        given:
        Global.session = Mock(Session) {
            getConfig() >> [tower: [
                endpoint: System.getenv('TOWER_ENDPOINT') ?: 'https://api.cloud.seqera.io',
                accessToken: System.getenv('TOWER_ACCESS_TOKEN'),
                workspaceId: System.getenv('TOWER_WORKSPACE_ID')
            ]]
        }

        and:
        def provider = new DatasetFileSystemProvider()
        def path = provider.getPath(new URI("dataset:///${datasetName}"))

        when:
        def bytes = provider.newInputStream(path).readNBytes(64)

        then:
        bytes.length > 0

        where:
        datasetName << datasetNames()
    }

    private static List<String> datasetNames() {
        System.getenv('DATASET_LIVE_NAMES')
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    }
}
