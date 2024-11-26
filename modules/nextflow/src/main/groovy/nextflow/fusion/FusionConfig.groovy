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

package nextflow.fusion


import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.util.MemoryUnit
/**
 * Model Fusion config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FusionConfig {

    final static public String DEFAULT_FUSION_AMD64_URL = 'https://fusionfs.seqera.io/releases/v2.4-amd64.json'
    final static public String DEFAULT_FUSION_ARM64_URL = 'https://fusionfs.seqera.io/releases/v2.4-arm64.json'
    final static public String DEFAULT_TAGS = "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"

    final static public String FUSION_PATH = '/usr/bin/fusion'

    final static private Pattern VERSION_JSON = ~/https:\/\/.*\/releases\/v(\d+(?:\.\w+)*)-(\w*)\.json$/

    final private Boolean enabled
    final private String containerConfigUrl
    @Deprecated final private Boolean exportAwsAccessKeys
    final private Boolean exportStorageCredentials
    final private String logOutput
    final private String logLevel
    final private boolean tagsEnabled
    final private String tagsPattern
    final private boolean privileged
    final private MemoryUnit cacheSize

    boolean enabled() { enabled }

    @Deprecated boolean exportAwsAccessKeys() { exportAwsAccessKeys }

    boolean exportStorageCredentials() {
        return exportStorageCredentials!=null
            ? exportStorageCredentials
            : exportAwsAccessKeys
    }

    String logLevel() { logLevel }

    String logOutput() { logOutput }

    boolean tagsEnabled() { tagsEnabled }

    String tagsPattern() { tagsPattern }

    MemoryUnit cacheSize() { cacheSize }

    URL containerConfigUrl() {
        this.containerConfigUrl ? new URL(this.containerConfigUrl) : null
    }

    boolean privileged() {
        return privileged
    }

    FusionConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled
        this.exportAwsAccessKeys = opts.exportAwsAccessKeys
        this.exportStorageCredentials = opts.exportStorageCredentials
        this.containerConfigUrl = opts.containerConfigUrl?.toString() ?: env.get('FUSION_CONTAINER_CONFIG_URL')
        this.logLevel = opts.logLevel
        this.logOutput = opts.logOutput
        this.tagsEnabled = opts.tags==null || opts.tags.toString()!='false'
        this.tagsPattern = (opts.tags==null || (opts.tags instanceof Boolean && opts.tags)) ? DEFAULT_TAGS : ( opts.tags !instanceof Boolean ? opts.tags as String : null )
        this.privileged = opts.privileged==null || opts.privileged.toString()=='true'
        this.cacheSize = opts.cacheSize as MemoryUnit
        if( containerConfigUrl && !validProtocol(containerConfigUrl))
            throw new IllegalArgumentException("Fusion container config URL should start with 'http:' or 'https:' protocol prefix - offending value: $containerConfigUrl")
    }

    protected boolean validProtocol(String url) {
        url.startsWith('http://') || url.startsWith('https://') || url.startsWith('file:/')
    }

    static FusionConfig getConfig() {
        return createConfig0(Global.config?.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
    }

    static FusionConfig getConfig(Session session) {
        return createConfig0(session.config?.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
    }

    @Memoized
    static private FusionConfig createConfig0(Map config, Map env) {
        new FusionConfig(config, env)
    }

    protected String retrieveFusionVersion(String url) {
        if( !url )
            return null
        final matcher_json = VERSION_JSON.matcher(url)
        if( matcher_json.matches() )
            return matcher_json.group(1)
        return null
    }

    String version() {
        return enabled
            ? retrieveFusionVersion(this.containerConfigUrl ?: DEFAULT_FUSION_AMD64_URL)
            : null
    }
}
