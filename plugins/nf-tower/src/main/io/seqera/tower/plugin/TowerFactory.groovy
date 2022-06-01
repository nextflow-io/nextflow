/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
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
        Boolean isEnabled = config.navigate('tower.enabled') as Boolean
        String endpoint = config.navigate('tower.endpoint') as String
        Duration requestInterval = config.navigate('tower.requestInterval') as Duration
        Duration aliveInterval = config.navigate('tower.aliveInterval') as Duration

        if( !isEnabled )
            return Collections.emptyList()

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
        final result = new ArrayList(1)
        result.add(tower)
        // register auth provider
        // note: this is needed to authorize access to resources via XFileSystemProvider used by NF
        // it's not needed by the tower client logic
        XAuthRegistry.instance.register(provider(tower.endpoint, tower.accessToken))
        return result
    }

    protected XAuthProvider provider(String endpoint, String accessToken) {
        assert !endpoint.endsWith('/'), "Tower endpoint URL should end with a `/` character"
        final pattern = ~/(?i)^$endpoint\/.*$/
        new XAuthProvider() {
            @Override
            boolean authorize(URLConnection conn) {
                final req = conn.getURL().toString()
                if( pattern.matcher(req).matches() && !conn.getRequestProperty('Authorization') ) {
                    log.trace "Authorizing request connection to: $req"
                    conn.setRequestProperty('Authorization', "Bearer $accessToken")
                    return true
                }
                return false
            }
        }
    }
}
