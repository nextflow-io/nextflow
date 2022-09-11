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
/**
 * Model Tower config accessed by Wave
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerConfig {

    final String accessToken

    final Long workspaceId

    TowerConfig(Map opts, Map<String,String> env) {
        this.accessToken = opts.accessToken as String ?: env.get('TOWER_ACCESS_TOKEN')
        this.workspaceId = opts.workspaceId as Long ?: env.get('TOWER_WORKSPACE_ID') as Long
    }
}
