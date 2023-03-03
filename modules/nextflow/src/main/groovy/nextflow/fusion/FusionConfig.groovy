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

package nextflow.fusion

import groovy.transform.CompileStatic
/**
 * Model Fusion config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FusionConfig {

    final static public String DEFAULT_FUSION_AMD64_URL = 'https://fusionfs.seqera.io/releases/v2.1-amd64.json'
    final static public String DEFAULT_FUSION_ARM64_URL = 'https://fusionfs.seqera.io/releases/v2.1-arm64.json'
    final static public String DEFAULT_TAGS = "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"

    final static public String FUSION_PATH = '/usr/bin/fusion'

    final private Boolean enabled
    final private String containerConfigUrl
    final private Boolean exportAwsAccessKeys
    final private String logOutput
    final private String logLevel
    final private boolean tagsEnabled
    final private String tagsPattern

    boolean enabled() { enabled }

    boolean exportAwsAccessKeys() { exportAwsAccessKeys }

    String logLevel() { logLevel }

    String logOutput() { logOutput }

    boolean tagsEnabled() { tagsEnabled }

    String tagsPattern() { tagsPattern }

    URL containerConfigUrl() {
        this.containerConfigUrl ? new URL(this.containerConfigUrl) : null
    }

    FusionConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled
        this.exportAwsAccessKeys = opts.exportAwsAccessKeys
        this.containerConfigUrl = opts.containerConfigUrl?.toString() ?: env.get('FUSION_CONTAINER_CONFIG_URL')
        this.logLevel = opts.logLevel
        this.logOutput = opts.logOutput
        this.tagsEnabled = opts.tags==null || opts.tags.toString()!='false'
        this.tagsPattern = (opts.tags==null || (opts.tags instanceof Boolean && opts.tags)) ? DEFAULT_TAGS : ( opts.tags !instanceof Boolean ? opts.tags as String : null )
        if( containerConfigUrl && !validProtocol(containerConfigUrl))
            throw new IllegalArgumentException("Fusion container config URL should start with 'http:' or 'https:' protocol prefix - offending value: $containerConfigUrl")
    }

    protected boolean validProtocol(String url) {
        url.startsWith('http://') || url.startsWith('https://') || url.startsWith('file:/')
    }

}
