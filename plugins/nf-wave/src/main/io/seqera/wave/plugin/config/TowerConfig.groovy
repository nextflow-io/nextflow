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

package io.seqera.wave.plugin.config

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Model Tower config accessed by Wave
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeNames = true, includePackage = false)
@CompileStatic
class TowerConfig {

    final String accessToken

    final String refreshToken

    final Long workspaceId

    final String endpoint

    final String workflowId

    TowerConfig(Map opts, Map<String,String> env) {
        this.accessToken = accessToken0(opts, env)
        this.refreshToken = refreshToken0(opts, env)
        this.workspaceId = workspaceId0(opts, env) as Long
        this.endpoint = endpoint0(opts, env)
        this.workflowId = env.get('TOWER_WORKFLOW_ID')
    }

    private String endpoint0(Map opts, Map<String,String> env) {
        def result = opts.endpoint as String
        if( !result || result=='-' )
            result = env.get('TOWER_API_ENDPOINT') ?: 'https://api.cloud.seqera.io'
        return result.stripEnd('/')
    }

    private String accessToken0(Map opts, Map<String,String> env) {
        // when 'TOWER_WORKFLOW_ID' is provided in the env, it's a tower made launch
        // therefore the access token should only be taken from the env
        // otherwise check into the config file and fallback in the env
        // see also
        // https://github.com/nextflow-io/nextflow/blob/master/plugins/nf-tower/src/main/io/seqera/tower/plugin/TowerClient.groovy#L369-L377
        def token = env.get('TOWER_WORKFLOW_ID')
                ? env.get('TOWER_ACCESS_TOKEN')
                : opts.containsKey('accessToken') ? opts.accessToken as String : env.get('TOWER_ACCESS_TOKEN')
        return token
    }

    private String refreshToken0(Map opts, Map<String,String> env) {
        // when 'TOWER_WORKFLOW_ID' is provided in the env, it's a tower made launch
        // therefore the access token should only be taken from the env
        // otherwise check into the config file and fallback in the env
        // see also
        // https://github.com/nextflow-io/nextflow/blob/master/plugins/nf-tower/src/main/io/seqera/tower/plugin/TowerClient.groovy#L369-L377
        def token = env.get('TOWER_WORKFLOW_ID')
                ? env.get('TOWER_REFRESH_TOKEN')
                : opts.containsKey('refreshToken') ? opts.refreshToken as String : env.get('TOWER_REFRESH_TOKEN')
        return token
    }

    private String workspaceId0(Map opts, Map<String,String> env) {
        // when 'TOWER_WORKFLOW_ID' is provided in the env, it's a tower made launch
        // therefore the workspace should only be taken from the env
        // otherwise check into the config file and fallback in the env
        def workspaceId = env.get('TOWER_WORKFLOW_ID')
                ? env.get('TOWER_WORKSPACE_ID')
                : opts.workspaceId as Long ?: env.get('TOWER_WORKSPACE_ID') as Long
        return workspaceId
    }
}
