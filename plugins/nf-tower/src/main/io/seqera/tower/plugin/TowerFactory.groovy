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

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.file.http.XAuthProvider
import nextflow.file.http.XAuthRegistry
import nextflow.platform.PlatformHelper
import nextflow.trace.TraceObserverFactoryV2
import nextflow.trace.TraceObserverV2
import nextflow.util.Duration
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
        final opts = session.config.tower as Map ?: Collections.emptyMap()
        final config = new TowerConfig(opts, env)
        if( !isEnabled(session, config, env) )
            return Collections.emptyList()
        // make sure the access token is available before the client is created, otherwise the
        // missing-token error is thrown deep in the TowerClient constructor as an AbortRunException
        // during session init and gets swallowed silently by the launcher
        checkAccessToken(config)
        final result = new ArrayList<TraceObserverV2>(1)
        final client = client(session, env)
        // resolve the workspace: a local setting wins, otherwise fall back to the
        // user's default workspace configured in Seqera Platform
        final local = PlatformHelper.getWorkspaceId(opts, env)
        final workspaceId = PlatformHelper.getEffectiveWorkspaceId(opts, env, () -> defaultWorkspaceId(opts))
        // when -- and only when -- the workspace came from the Platform account default,
        // publish it back into the session config. Other subsystems that scope themselves
        // to the workspace -- Wave (registry credentials), the Fusion licence, the Seqera
        // executor -- read it via PlatformHelper.getWorkspaceId(session.config.tower, env)
        // and would otherwise resolve to the personal workspace while the run is reported
        // into the Platform default one. Anything the user set is left exactly as it is.
        if( !local && workspaceId )
            publishDefaultWorkspaceId(session, workspaceId)
        // create the tower observer
        result.add( new TowerObserver(session, client, workspaceId, env))
        // create the logs checkpoint
        if( session.cloudCachePath )
            result.add( new LogsCheckpoint() )
        return result
    }

    /**
     * Query the user's server-side default workspace from Seqera Platform.
     *
     * This runs during session initialization, before the pipeline script is even parsed,
     * so it uses a dedicated client with a bounded retry budget and short timeouts: the
     * lookup is a best-effort convenience that falls back to the personal workspace, and
     * it must never be able to stall the start of a run for minutes when Platform is slow
     * or unreachable. The shared client's retry policy is sized for the telemetry stream,
     * which is a different trade-off.
     *
     * Extracted as a seam so it can be stubbed in tests without hitting the network.
     */
    protected String defaultWorkspaceId(Map opts) {
        return new TowerClient(TowerConfig.forLookup(opts, env)).getDefaultWorkspaceId()
    }

    /**
     * Store the workspace resolved from the Seqera Platform account default in the session
     * config, so that every subsystem scoping itself to the workspace observes it.
     *
     * Only ever called when the user set no workspace of their own, so this never
     * overwrites a configured value.
     *
     * Note this relies on {@code session.config} being a plain map: it is normalized from
     * the parsed {@code ConfigObject} by {@code ConfigCmdAdapter.normalize0} and converted
     * in the {@link Session} constructor, so an absent key really does read as null here.
     */
    protected void publishDefaultWorkspaceId(Session session, String workspaceId) {
        final config = session.config
        if( config.tower == null )
            config.tower = new HashMap(1)
        (config.tower as Map).workspaceId = workspaceId
        // tell the user why their run is landing in a workspace they did not configure
        log.info "Using default workspace configured in your Seqera Platform account: $workspaceId"
    }

    @Memoized
    static TowerClient client(Session session, Map env) {
        final opts = session.config.tower as Map ?: Collections.emptyMap()
        final config = new TowerConfig(opts, env)
        final tower = new TowerClient(config)
        // register auth provider
        // note: this is needed to authorize access to resources via XFileSystemProvider used by NF
        // it's not needed by the tower client logic
        XAuthRegistry.instance.register(provider(config.endpoint, config.accessToken))
        return tower
    }

    static protected XAuthProvider provider(String endpoint, String accessToken) {
        if (endpoint.endsWith('/'))
            throw new IllegalArgumentException("Seqera Platform endpoint URL should not end with a `/` character -- offending value: $endpoint")
        final refreshToken = SysEnv.get('TOWER_REFRESH_TOKEN')
        return new TowerXAuth(endpoint, accessToken, refreshToken)
    }

    private static boolean isEnabled(Session session, TowerConfig config, Map<String,String> env) {
        return config.enabled || env.get('TOWER_WORKFLOW_ID') || session.config.navigate('fusion.enabled') as Boolean
    }

    /**
     * Verify a Seqera Platform access token is available before creating the observer.
     *
     * The Tower observer can be enabled solely because Fusion is enabled (see {@link #isEnabled}).
     * Checking the token here and raising an {@link AbortOperationException} reports an actionable
     * message to the user, instead of the silent {@code AbortRunException} thrown later by the
     * {@link TowerClient} constructor.
     */
    protected void checkAccessToken(TowerConfig config) {
        if( config.accessToken )
            return
        throw new AbortOperationException("Seqera Platform access token is required -- Provide it via the `tower.accessToken` config setting or the `TOWER_ACCESS_TOKEN` environment variable")
    }

    static TowerClient client() {
        return client(Global.session as Session, SysEnv.get())
    }
}
