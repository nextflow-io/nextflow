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

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.file.http.XAuthProvider
import nextflow.file.http.XAuthRegistry
import nextflow.trace.TraceObserverV2
import nextflow.trace.TraceObserverFactoryV2
import nextflow.util.Duration
import nextflow.util.SimpleHttpClient
/**
 * Create and register the Tower observer instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerFactory implements TraceObserverFactoryV2 {

    private Map<String,String> env

    TowerFactory(){
        env = SysEnv.get()
    }

    @Override
    Collection<TraceObserverV2> create(Session session) {
        final client = client(session, env)
        if( !client )
            return Collections.emptyList()

        final result = new ArrayList<TraceObserverV2>(1)
        // create the tower client
        result.add(client)
        // create the logs checkpoint
        if( session.cloudCachePath )
            result.add( new LogsCheckpoint() )
        return result
    }

    static protected TowerClient createTowerClient0(Session session, Map config, Map env) {
        String endpoint = config.navigate('tower.endpoint') as String
        Duration requestInterval = config.navigate('tower.requestInterval') as Duration
        Duration aliveInterval = config.navigate('tower.aliveInterval') as Duration

        if ( !endpoint || endpoint=='-' )
            endpoint = env.get('TOWER_API_ENDPOINT') ?: TowerClient.DEF_ENDPOINT_URL

        final tower = new TowerClient(session, endpoint).withEnvironment(env)
        if( aliveInterval )
            tower.aliveInterval = aliveInterval
        if( requestInterval )
            tower.requestInterval = requestInterval
        // error handling settings
        tower.maxRetries = config.navigate('tower.maxRetries', 5) as int
        tower.backOffBase = config.navigate('tower.backOffBase', SimpleHttpClient.DEFAULT_BACK_OFF_BASE) as int
        tower.backOffDelay = config.navigate('tower.backOffDelay', SimpleHttpClient.DEFAULT_BACK_OFF_DELAY  ) as int
        // when 'TOWER_WORKFLOW_ID' is provided in the env, it's a tower made launch
        // therefore the workspace should only be taken from the env
        // otherwise check into the config file and fallback in the env
        tower.workspaceId = env.get('TOWER_WORKFLOW_ID')
                ? env.get('TOWER_WORKSPACE_ID')
                : config.navigate('tower.workspaceId', env.get('TOWER_WORKSPACE_ID'))

        // register auth provider
        // note: this is needed to authorize access to resources via XFileSystemProvider used by NF
        // it's not needed by the tower client logic
        XAuthRegistry.instance.register(provider(tower.endpoint, tower.accessToken))
        return tower
    }

    static protected XAuthProvider provider(String endpoint, String accessToken) {
        if (endpoint.endsWith('/'))
            throw new IllegalArgumentException("Seqera Platform endpoint URL should not end with a `/` character -- offending value: $endpoint")
        final refreshToken = SysEnv.get('TOWER_REFRESH_TOKEN')
        return new TowerXAuth(endpoint, accessToken, refreshToken)
    }

    @Memoized
    static TowerClient client(Session session, Map<String,String> env) {
        final config = session.config
        Boolean isEnabled = config.navigate('tower.enabled') as Boolean || env.get('TOWER_WORKFLOW_ID') || config.navigate('fusion.enabled') as Boolean
        return isEnabled
            ? createTowerClient0(session, config, env)
            : null
    }

    static TowerClient client() {
        client(Global.session as Session, SysEnv.get())
    }
}
