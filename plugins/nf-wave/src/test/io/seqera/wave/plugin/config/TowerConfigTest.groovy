/*
 * Copyright 2020-2022, Seqera Labs
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

package io.seqera.wave.plugin.config

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerConfigTest extends Specification {

    def 'should get tower config' () {
        given:
        TowerConfig config
        Map env

        when:
        config = new TowerConfig([:], [:])
        then:
        !config.getAccessToken()
        !config.getWorkspaceId()

        when:
        env = [TOWER_ACCESS_TOKEN:'foo', TOWER_WORKSPACE_ID: '123']
        config = new TowerConfig([:], env)
        then:
        config.accessToken == 'foo'
        config.workspaceId == 123

        when:
        env =  [TOWER_ACCESS_TOKEN:'foo', TOWER_WORKSPACE_ID: '123']
        config = new TowerConfig(accessToken: 'bar', workspaceId: '789', env)
        then:
        config.accessToken == 'bar'
        config.workspaceId == 789

        when:
        env =  [TOWER_ACCESS_TOKEN:'foo', TOWER_WORKSPACE_ID: '123']
        config = new TowerConfig(accessToken: null, workspaceId: '789', env)
        then:
        config.accessToken == null
        config.workspaceId == 789

        // when TOWER_WORKFLOW_ID is defined env has priority
        when:
        env = [TOWER_ACCESS_TOKEN:'foo', TOWER_WORKSPACE_ID: '123', TOWER_WORKFLOW_ID: 'xyz']
        config = new TowerConfig(accessToken: 'bar', workspaceId: '789', env)
        then:
        config.accessToken == 'foo'
        config.workspaceId == 123
    }
}
