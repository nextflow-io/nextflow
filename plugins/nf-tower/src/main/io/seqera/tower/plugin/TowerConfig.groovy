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
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.ScopeName
import nextflow.script.dsl.Description
import nextflow.platform.PlatformHelper

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

    @ConfigOption
    @Description("""
        The unique access token for your Seqera Platform account.
    """)
    final String accessToken

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
        The workspace ID in Seqera Platform in which to save the run (default: the launching user's personal workspace).
    """)
    final String workspaceId

    final TowerRetryPolicy retryPolicy

    /* required by extension point -- do not remove */
    TowerConfig() {}

    TowerConfig(Map opts, Map<String,String> env) {
        this.accessToken = PlatformHelper.getAccessToken(opts, env)
        this.enabled = opts.enabled as boolean
        this.endpoint = PlatformHelper.getEndpoint(opts, env)
        this.workspaceId = PlatformHelper.getWorkspaceId(opts, env)
        this.retryPolicy = new TowerRetryPolicy(opts.retryPolicy as Map ?: Map.of(), opts)
    }
}
