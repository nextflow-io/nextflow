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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.http.XAuthProvider
import nextflow.file.http.XAuthRegistry
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory
import nextflow.util.Duration
import nextflow.util.SimpleHttpClient
/**
 * Create and register the Tower observer instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerFactory implements TraceObserverFactory {

    private Map<String,String> env

    TowerFactory(){
        env = System.getenv()
    }

    @Override
    Collection<TraceObserver> create(Session session) {
        final config = session.config
        Boolean isEnabled = config.navigate('tower.enabled') as Boolean || env.get('TOWER_WORKFLOW_ID')

        if( !isEnabled )
            return Collections.emptyList()

        final result = new ArrayList(1)
        // create the tower client
        final tower = createTowerClient(session, config)
        result.add(tower)
        // create the logs checkpoint
        if( session.cloudCachePath )
            result.add( new LogsCheckpoint() )
        return result
    }

    protected TowerClient createTowerClient(Session session, Map config) {
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

    protected XAuthProvider provider(String endpoint, String accessToken) {
        if (endpoint.endsWith('/'))
            throw new IllegalArgumentException("Seqera Platform endpoint URL should not end with a `/` character -- offending value: $endpoint")
        final refreshToken = env.get('TOWER_REFRESH_TOKEN')
        return new TowerXAuth(endpoint, accessToken, refreshToken)
    }

}
