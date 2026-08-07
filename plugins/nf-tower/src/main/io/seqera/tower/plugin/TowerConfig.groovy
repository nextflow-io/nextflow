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
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.platform.PlatformHelper
import nextflow.util.Duration

/**
 * Model Seqera Platform configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("tower")
@Description("""
    The `tower` scope controls the settings for [Seqera Platform](https://seqera.io) (formerly Tower Cloud).
""")
@CompileStatic
class TowerConfig implements ConfigScope {

    static final Duration DEFAULT_CONNECT_TIMEOUT = Duration.of('60s')

    static final Duration DEFAULT_READ_TIMEOUT = Duration.of('60s')

    /**
     * Timeout applied to a best-effort lookup. Note this is a *per-phase* budget: it is
     * used for both the connect and the read phase, so the worst case for one request is
     * roughly twice this value.
     */
    static final Duration LOOKUP_PHASE_TIMEOUT = Duration.of('10s')

    @ConfigOption
    @Description("""
        The unique access token for your Seqera Platform account.
    """)
    final String accessToken

    @ConfigOption
    @Description("""
        The Compute Environment ID in Seqera Platform in which to launch the run (default: the primary environment in the workspace).
    """)
    final String computeEnvId

    @ConfigOption
    @Description("""
        Enable workflow monitoring with Seqera Platform (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The endpoint of your Seqera Platform instance (default: `https://api.cloud.seqera.io`).
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        The HTTP connection timeout for Seqera Platform API requests (default: `'60s'`).
    """)
    final Duration httpConnectTimeout

    @ConfigOption
    @Description("""
        The HTTP read timeout for Seqera Platform API requests (default: `'60s'`).
    """)
    final Duration httpReadTimeout

    final TowerRetryPolicy retryPolicy

    @ConfigOption
    @Description("""
        The workspace ID in Seqera Platform in which to save the run (default: the user's default workspace in Seqera Platform, or the launching user's personal workspace if no default is set).
    """)
    final String workspaceId

    /**
     * Build a config for a best-effort lookup: one that improves the experience when it
     * succeeds but must never delay or abort the operation when Platform is slow or
     * unreachable, because the caller has a working fallback.
     *
     * The default retry policy is sized for the telemetry stream -- 10 attempts with
     * exponential backoff, up to roughly three minutes -- which is the wrong trade-off
     * for such a lookup, especially on the session-init path where it delays the start of
     * every run. This bounds it to a single attempt with short timeouts.
     *
     * The caller's options are copied, so the map passed in is never modified.
     *
     * @param opts the `tower` scope options, keyed without the prefix
     * @param env the applicable environment variables
     * @return a {@link TowerConfig} identical to the caller's except for retries and timeouts
     */
    static TowerConfig forLookup(Map opts, Map<String,String> env) {
        final bounded = new HashMap(opts)
        bounded.retryPolicy = [maxAttempts: 1]
        bounded.httpConnectTimeout = LOOKUP_PHASE_TIMEOUT
        bounded.httpReadTimeout = LOOKUP_PHASE_TIMEOUT
        return new TowerConfig(bounded, env)
    }

    /* required by extension point -- do not remove */
    TowerConfig() {}

    TowerConfig(Map opts, Map<String,String> env) {
        this.accessToken = PlatformHelper.getAccessToken(opts, env)
        if( opts.computeEnvId )
            this.computeEnvId = opts.computeEnvId as String
        this.enabled = opts.enabled as boolean
        this.endpoint = PlatformHelper.getEndpoint(opts, env)
        this.httpConnectTimeout = opts.httpConnectTimeout ? opts.httpConnectTimeout as Duration : DEFAULT_CONNECT_TIMEOUT
        this.httpReadTimeout = opts.httpReadTimeout ? opts.httpReadTimeout as Duration : DEFAULT_READ_TIMEOUT
        this.retryPolicy = new TowerRetryPolicy(opts.retryPolicy as Map ?: Map.of(), opts)
        this.workspaceId = PlatformHelper.getWorkspaceId(opts, env)
    }
}
