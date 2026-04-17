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

package io.seqera.executor

import io.seqera.config.SeqeraConfig
import io.seqera.sched.client.SchedClientConfig
import nextflow.Session
import nextflow.SysEnv
import nextflow.platform.PlatformHelper
import spock.lang.Specification

/**
 * Tests for SeqeraExecutor client configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SeqeraExecutorTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'should create client config with config settings'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [endpoint: 'https://api.platform.example.com', accessToken: 'config-access-token', refreshToken: 'config-refresh-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.platform.example.com'
        config.accessToken == 'config-access-token'
        config.refreshToken == 'config-refresh-token'
    }

    def 'should create client config with env variable settings'() {
        given:
        SysEnv.push([
            TOWER_API_ENDPOINT: 'https://api.env.example.com',
            TOWER_ACCESS_TOKEN: 'env-access-token',
            TOWER_REFRESH_TOKEN: 'env-refresh-token'
        ])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [:]
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.env.example.com'
        config.accessToken == 'env-access-token'
        config.refreshToken == 'env-refresh-token'
    }

    def 'should use default platform url when not configured'() {
        given:
        SysEnv.push([:])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [accessToken: 'my-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.cloud.seqera.io'
        config.accessToken == 'my-token'
        config.refreshToken == null
    }

    def 'should prefer config over env variables'() {
        given:
        SysEnv.push([
            TOWER_API_ENDPOINT: 'https://api.env.example.com',
            TOWER_ACCESS_TOKEN: 'env-access-token',
            TOWER_REFRESH_TOKEN: 'env-refresh-token'
        ])

        when:
        def config = buildClientConfig(
            [endpoint: 'https://sched.example.com', region: 'us-west-2'],
            [endpoint: 'https://api.config.example.com', accessToken: 'config-access-token', refreshToken: 'config-refresh-token']
        )

        then:
        config.endpoint == 'https://sched.example.com'
        config.platformUrl == 'https://api.config.example.com'
        config.accessToken == 'config-access-token'
        config.refreshToken == 'config-refresh-token'
    }

    def 'should set fusion default version when not configured' () {
        given:
        SysEnv.push([:])
        def fusionConfig = [enabled: true]
        def config = [fusion: fusionConfig]
        def session = Mock(Session) { getConfig() >> config }
        def executor = new SeqeraExecutor(session: session)

        when:
        executor.applyFusionDefaults()

        then:
        fusionConfig.targetVersion == '2.6'
    }

    def 'should not override fusion version when containerConfigUrl is set' () {
        given:
        SysEnv.push([:])
        def fusionConfig = [enabled: true, containerConfigUrl: 'https://custom.url/v3.0-amd64.json']
        def config = [fusion: fusionConfig]
        def session = Mock(Session) { getConfig() >> config }
        def executor = new SeqeraExecutor(session: session)

        when:
        executor.applyFusionDefaults()

        then:
        fusionConfig.targetVersion == null
    }

    def 'should expose run resource labels coerced from config-level process.resourceLabels'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        executor.session = Mock(Session) {
            getConfig() >> [process: [resourceLabels: [team: 'a', priority: 7]]]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == [team: 'a', priority: '7']
    }

    def 'should yield empty run resource labels when process.resourceLabels is absent'() {
        given:
        SysEnv.push([:])
        def executor = new SeqeraExecutor()
        executor.session = Mock(Session) {
            getConfig() >> [:]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == [:]
    }

    /**
     * Builds a SchedClientConfig using the same logic as {@link SeqeraExecutor#createClient()}
     */
    private SchedClientConfig buildClientConfig(Map executorOpts, Map towerConfig) {
        def seqeraConfig = new SeqeraConfig([executor: executorOpts]).executor
        def accessToken = PlatformHelper.getAccessToken(towerConfig, SysEnv.get())
        def refreshToken = PlatformHelper.getRefreshToken(towerConfig, SysEnv.get())
        def platformUrl = PlatformHelper.getEndpoint(towerConfig, SysEnv.get())
        return SchedClientConfig.builder()
                .endpoint(seqeraConfig.endpoint)
                .platformUrl(platformUrl)
                .accessToken(accessToken)
                .refreshToken(refreshToken)
                .retryConfig(seqeraConfig.retryOpts())
                .build()
    }
}
