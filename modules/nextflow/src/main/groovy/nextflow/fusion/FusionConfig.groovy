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

package nextflow.fusion


import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.util.MemoryUnit
/**
 * Model Fusion config options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("fusion")
@Description("""
    The `fusion` scope provides advanced configuration for the use of the [Fusion file system](https://docs.seqera.io/fusion).
""")
@CompileStatic
class FusionConfig implements ConfigScope {

    final static public String DEFAULT_FUSION_AMD64_URL = 'https://fusionfs.seqera.io/releases/v2.5-amd64.json'
    final static public String DEFAULT_FUSION_ARM64_URL = 'https://fusionfs.seqera.io/releases/v2.5-arm64.json'
    final static public String DEFAULT_SNAPSHOT_AMD64_URL = 'https://fusionfs.seqera.io/releases/v2.5-snap_amd64.json'
    final static public String DEFAULT_SNAPSHOT_ARM64_URL = 'https://fusionfs.seqera.io/releases/v2.5-snap_arm64.json'

    final static public String DEFAULT_TAGS = "[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)"

    final static public int DEFAULT_SNAPSHOT_MAX_SPOT_ATTEMPTS = 5

    final static public String FUSION_PATH = '/usr/bin/fusion'

    final static private String PRODUCT_NAME = 'fusion'

    final static private Pattern VERSION_JSON = ~/https:\/\/.*\/releases\/v(\d+(?:\.\w+)*)-(\w*)\.json$/

    @ConfigOption
    @Description("""
        Enable the Fusion file system (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The maximum size of the local cache used by the Fusion client.
    """)
    final MemoryUnit cacheSize

    @ConfigOption
    @Description("""
        The URL of the container layer that provides the Fusion client.
    """)
    final String containerConfigUrl

    @ConfigOption
    @Description("""
        Export the access credentials required by the underlying object storage to the task execution environment (default: `false`).
    """)
    final boolean exportStorageCredentials

    @ConfigOption
    @Description("""
        The log level of the Fusion client.
    """)
    final String logLevel

    @ConfigOption
    @Description("""
        The output location of the Fusion log.
    """)
    final String logOutput

    @ConfigOption
    @Description("""
        Enable privileged containers for Fusion (default: `true`).
    """)
    final boolean privileged

    @ConfigOption
    @Description("""
        Enable Fusion snapshotting (preview, default: `false`). This feature allows Fusion to automatically restore a job when it is interrupted by a spot reclamation.
    """)
    final boolean snapshots

    @ConfigOption(types=[Boolean])
    @Description("""
        The pattern that determines how tags are applied to files created via the Fusion client (default: `[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true)`). Set to `false` to disable tags.
    """)
    final String tags

    boolean enabled() { enabled }

    boolean exportStorageCredentials() {
        return exportStorageCredentials
    }

    String logLevel() { logLevel }

    String logOutput() { logOutput }

    boolean tagsEnabled() { tags != null }

    String tagsPattern() { tags }

    MemoryUnit cacheSize() { cacheSize }

    boolean snapshotsEnabled() { snapshots }

    URL containerConfigUrl() {
        this.containerConfigUrl ? new URL(this.containerConfigUrl) : null
    }

    boolean privileged() {
        return privileged
    }

    /* required by extension point -- do not remove */
    FusionConfig() {}

    FusionConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled as boolean
        this.exportStorageCredentials = (opts.exportStorageCredentials ?: opts.exportAwsAccessKeys) as boolean
        this.containerConfigUrl = opts.containerConfigUrl ?: env.get('FUSION_CONTAINER_CONFIG_URL')
        this.logLevel = opts.logLevel
        this.logOutput = opts.logOutput
        this.tags = parseTags(opts.tags)
        this.privileged = opts.privileged == null || opts.privileged as boolean
        this.cacheSize = opts.cacheSize as MemoryUnit
        this.snapshots = opts.snapshots as boolean

        if( containerConfigUrl && !validProtocol(containerConfigUrl))
            throw new IllegalArgumentException("Fusion container config URL should start with 'http:' or 'https:' protocol prefix - offending value: $containerConfigUrl")
    }

    static private String parseTags(Object value) {
        if( value == null )
            return DEFAULT_TAGS
        if( value instanceof Boolean && value )
            return DEFAULT_TAGS
        if( value instanceof CharSequence )
            return value
        return null
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

    /**
     * Return the Fusion SKU string
     *
     * @return A string representing the Fusion SKU
     */
    String sku() {
        return enabled ? PRODUCT_NAME : null
    }

    String version() {
        return enabled
            ? retrieveFusionVersion(this.containerConfigUrl ?: DEFAULT_FUSION_AMD64_URL)
            : null
    }
}
