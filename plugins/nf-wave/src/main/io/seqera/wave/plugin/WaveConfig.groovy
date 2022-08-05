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

package io.seqera.wave.plugin

import groovy.transform.CompileStatic
import nextflow.util.Duration

/**
 * Model Wave client configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WaveConfig {
    final private static String DEF_ENDPOINT = 'http://localhost:9090'
    final private Boolean enabled
    final private String endpoint
    final private String containerConfigUrl
    final private Duration tokensCacheMaxDuration

    WaveConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled
        this.endpoint = (opts.endpoint?.toString() ?: env.get('WAVE_API_ENDPOINT') ?: DEF_ENDPOINT)?.stripEnd('/')
        this.containerConfigUrl = (opts.containerConfigUrl?.toString() ?: env.get('WAVE_CONTAINER_CONFIG_URL'))?.stripEnd('/')
        this.tokensCacheMaxDuration = opts.navigate('tokens.cache.maxDuration', '15m') as Duration
    }

    Boolean enabled() { this.enabled }

    String endpoint() { this.endpoint }

    Duration tokensCacheMaxDuration() { return tokensCacheMaxDuration }

    URL containerConfigUrl() {
        this.containerConfigUrl ? new URL(this.containerConfigUrl) : null
    }
}
