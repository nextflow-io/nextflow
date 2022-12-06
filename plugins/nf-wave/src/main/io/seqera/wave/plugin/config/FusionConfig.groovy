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

import groovy.transform.CompileStatic

/**
 * Model Fusion config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FusionConfig {

    final static public String DEFAULT_FUSION_AMD64_URL = 'https://fusionfs.seqera.io/releases/v0.6.0-amd64.json'
    final static public String DEFAULT_FUSION_ARM64_URL = 'https://fusionfs.seqera.io/releases/v0.6.0-arm64.json'

    final private enabled
    final private String containerConfigUrl

    boolean enabled() { enabled }

    URL containerConfigUrl() {
        this.containerConfigUrl ? new URL(this.containerConfigUrl) : null
    }

    FusionConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled
        this.containerConfigUrl = opts.containerConfigUrl?.toString() ?: env.get('FUSION_CONTAINER_CONFIG_URL')
        if( containerConfigUrl && (!containerConfigUrl.startsWith('http://') && !containerConfigUrl.startsWith('https://')))
            throw new IllegalArgumentException("Fusion container config URL should start with 'http:' or 'https:' protocol prefix - offending value: $containerConfigUrl")
    }

}
