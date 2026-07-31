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

package io.seqera.wave.plugin.config

import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import io.seqera.wave.api.BuildCompression
import io.seqera.wave.api.ScanLevel
import io.seqera.wave.api.ScanMode
import io.seqera.wave.config.CondaOpts
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description
import nextflow.file.FileHelper
import nextflow.util.Duration
/**
 * Model Wave client configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("wave")
@Description("""
    The `wave` scope provides advanced configuration for the use of [Wave containers](https://docs.seqera.io/wave).
""")
@Slf4j
@ToString(includeNames = true, includePackage = false, includeFields = true, useGetters = false)
@CompileStatic
class WaveConfig implements ConfigScope {

    final private static String DEF_ENDPOINT = 'https://wave.seqera.io'

    final private static List<String> DEF_STRATEGIES = List.of('container','dockerfile','conda')

    final WaveBuildConfig build

    @ConfigOption
    @Description("""
        Enable the use of Wave containers (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        The Wave service endpoint (default: `https://wave.seqera.io`).
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        Enable Wave container freezing (default: `false`). Wave will provision a non-ephemeral container image that will be pushed to a container repository of your choice.

        See also: `wave.build.repository` and `wave.build.cacheRepository`
    """)
    final boolean freeze

    final HttpOpts httpClient

    @ConfigOption
    @Description("""
        Enable Wave container mirroring (default: `false`). Wave will mirror (i.e. copy) the containers in your pipeline to a container registry of your choice, so that pipeline tasks can pull the containers from this registry instead of the original one.

        See also: `wave.build.repository`
    """)
    final boolean mirror

    final RetryOpts retryPolicy

    final WaveScanConfig scan

    @ConfigOption(types=[String])
    @Description("""
        The strategy to be used when resolving multiple Wave container requirements (default: `'container,dockerfile,conda'`).
    """)
    final List<String> strategy

    final private Boolean bundleProjectResources
    final private List<URL> containerConfigUrl
    final private Boolean preserveFileTimestamp
    final private Duration tokensCacheMaxDuration

    /* required by extension point -- do not remove */
    WaveConfig() {}

    WaveConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.build = new WaveBuildConfig(opts.build as Map ?: Collections.emptyMap())
        this.enabled = opts.enabled as boolean
        this.endpoint = (opts.endpoint?.toString() ?: env.get('WAVE_API_ENDPOINT') ?: DEF_ENDPOINT)?.stripEnd('/')
        this.freeze = opts.freeze as boolean
        this.httpClient = new HttpOpts(opts.httpClient as Map ?: Collections.emptyMap())
        this.mirror = opts.mirror as boolean
        this.retryPolicy = retryOpts0(opts)
        this.scan = new WaveScanConfig(opts.scan as Map ?: Collections.emptyMap())
        this.strategy = parseStrategy(opts.strategy)

        this.bundleProjectResources = opts.bundleProjectResources
        this.containerConfigUrl = parseConfig(opts, env)
        this.preserveFileTimestamp = opts.preserveFileTimestamp as Boolean
        this.tokensCacheMaxDuration = opts.navigate('tokens.cache.maxDuration', '30m') as Duration

        validateConfig()
    }

    Boolean enabled() { this.enabled }

    String endpoint() { this.endpoint }

    CondaOpts condaOpts() { this.build.conda.toCondaOpts() }

    RetryOpts retryOpts() { this.retryPolicy }

    HttpOpts httpOpts() { this.httpClient }

    List<String> strategy() { this.strategy }

    boolean freezeMode() { this.freeze }

    boolean mirrorMode() { this.mirror }

    boolean preserveFileTimestamp() { return this.preserveFileTimestamp }

    boolean bundleProjectResources() { bundleProjectResources }

    String buildRepository() { build.repository }

    String cacheRepository() { build.cacheRepository }

    String buildTemplate() { build.template }

    Duration buildMaxDuration() { build.maxDuration }

    BuildCompression buildCompression() { build.compression?.toBuildCompression() }

    private void validateConfig() {
        def scheme= FileHelper.getUrlProtocol(endpoint)
        if( scheme !in ['http','https'] )
            throw new IllegalArgumentException("Endpoint URL should start with 'http:' or 'https:' protocol prefix - offending value: '$endpoint'")
        if( FileHelper.getUrlProtocol(build.repository) )
            throw new IllegalArgumentException("Config setting 'wave.build.repository' should not include any protocol prefix - offending value: '$build.repository'")
        if( FileHelper.getUrlProtocol(build.cacheRepository) )
            throw new IllegalArgumentException("Config setting 'wave.build.cacheRepository' should not include any protocol prefix - offending value: '$build.cacheRepository'")
    }

    private RetryOpts retryOpts0(Map opts) {
        if( opts.retryPolicy )
            return new RetryOpts(opts.retryPolicy as Map)
        if( opts.retry ) {
            log.warn "Configuration options 'wave.retry' has been deprecated - replace it with 'wave.retryPolicy'"
            return new RetryOpts(opts.retry as Map)
        }
        return new RetryOpts(Collections.emptyMap())
    }

    protected List<String> parseStrategy(value) {
        if( !value ) {
            log.debug "Wave strategy not specified - using default: $DEF_STRATEGIES"
            return DEF_STRATEGIES
        }
        List<String> result
        if( value instanceof CharSequence )
            result = value.tokenize(',') .collect(it -> it.toString().trim())
        else if( value instanceof List )
            result = value.collect(it -> it.toString().trim())
        else
            throw new IllegalArgumentException("Invalid value for 'wave.strategy' configuration attribute - offending value: $value")
        for( String it : result ) {
            if( it !in DEF_STRATEGIES)
                throw new IllegalArgumentException("Invalid value for 'wave.strategy' configuration attribute - offending value: $it")
        }
        return result
    }

    protected List<URL> parseConfig(Map opts, Map<String,String> env) {
        List<String> result = new ArrayList<>(10)
        if( !opts.containerConfigUrl && env.get('WAVE_CONTAINER_CONFIG_URL') ) {
            result.add(checkUrl(env.get('WAVE_CONTAINER_CONFIG_URL')))
        }
        else if( opts.containerConfigUrl instanceof CharSequence ) {
            result.add(checkUrl(opts.containerConfigUrl.toString()))
        }
        else if( opts.containerConfigUrl instanceof List ) {
            for( def it : opts.containerConfigUrl ) {
                result.add(checkUrl(it.toString()))
            }
        }

        return result.collect(it -> new URL(it))
    }

    private String checkUrl(String value) {
        if( value && (!value.startsWith('http://') && !value.startsWith('https://')))
            throw new IllegalArgumentException("Wave container config URL should start with 'http:' or 'https:' protocol prefix - offending value: $value")
        return value
    }

    List<URL> containerConfigUrl() {
        return containerConfigUrl ?: Collections.<URL>emptyList()
    }

    Duration tokensCacheMaxDuration() {
        return tokensCacheMaxDuration
    }

    ScanMode scanMode() {
        return scan.mode
    }

    List<ScanLevel> scanAllowedLevels() {
        return scan.allowedLevels
    }
}


@ToString(includeNames = true, includePackage = false, includeFields = true)
@CompileStatic
class WaveScanConfig implements ConfigScope {

    @ConfigOption(types=[String])
    @Description("""
        Comma-separated list of allowed vulnerability levels when scanning containers for security vulnerabilities in `required` mode.
    """)
    final List<ScanLevel> allowedLevels

    @ConfigOption(types=[String])
    @Description("""
        Enable Wave container security scanning. Wave will scan the containers in your pipeline for security vulnerabilities.
    """)
    final ScanMode mode

    WaveScanConfig(Map opts) {
        allowedLevels = parseScanLevels(opts.allowedLevels)
        mode = opts.mode as ScanMode
    }

    protected List<ScanLevel> parseScanLevels(value) {
        if( !value )
            return null
        if( value instanceof CharSequence ) {
            final str = value.toString()
            value = str.tokenize(',').collect(it->it.trim())
        }
        if( value instanceof List ) {
            return (value as List).collect(it-> ScanLevel.valueOf(it.toString().toUpperCase()))
        }
        throw new IllegalArgumentException("Invalid value for 'wave.scan.levels' setting - offending value: $value; type: ${value.getClass().getName()}")
    }
}


@ToString(includeNames = true, includePackage = false, includeFields = true)
@CompileStatic
class WaveBuildConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The container repository where images built by Wave are uploaded.
    """)
    final String repository

    @ConfigOption
    @Description("""
        The container repository used to cache image layers built by the Wave service.
    """)
    final String cacheRepository

    @ConfigOption
    @Description("""
        The build template to use for container builds. Supported values: `conda/pixi:v1` (Pixi with multi-stage builds), `conda/micromamba:v2` (Micromamba 2.x with multi-stage builds), `cran/installr:v1` (R/CRAN packages). Default: standard conda/micromamba:v1 template.
    """)
    final String template

    @Description("""
        The `wave.build.conda` scope controls how Conda packages are provisioned into containers.
    """)
    final WaveBuildCondaConfig conda

    @Description("""
        The `wave.build.compression` scope controls compression of the container image built by Wave.
    """)
    final WaveBuildCompressionConfig compression

    @ConfigOption
    @Description("""
    """)
    final Duration maxDuration

    WaveBuildConfig(Map opts) {
        repository = opts.repository
        cacheRepository = opts.cacheRepository
        template = opts.template
        conda = new WaveBuildCondaConfig(opts.conda as Map ?: Collections.emptyMap())
        compression = new WaveBuildCompressionConfig(opts.compression as Map ?: Collections.emptyMap())
        maxDuration = opts.maxDuration as Duration ?: Duration.of('40m')
    }
}


/**
 * Model the `wave.build.conda` config scope.
 *
 * This is a Nextflow-side mirror of the wave-api {@link CondaOpts} type, which
 * cannot implement {@link ConfigScope} because it belongs to an external
 * library. Use {@link #toCondaOpts()} to obtain the corresponding domain object.
 */
@ToString(includeNames = true, includePackage = false, includeFields = true)
@CompileStatic
class WaveBuildCondaConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The Mamba container image used to build Conda based containers. This is expected to be a [micromamba-docker](https://github.com/mamba-org/micromamba-docker) image.
    """)
    final String mambaImage

    @ConfigOption
    @Description("""
        The base image for the final stage in multi-stage Conda container builds (default: `ubuntu:24.04`). This option only applies when `wave.build.template` is set to `conda/micromamba:v2` or `conda/pixi:v1`.
    """)
    final String baseImage

    @ConfigOption
    @Description("""
        One or more Conda packages to be always added in the resulting container (default: `conda-forge::procps-ng`).
    """)
    final String basePackages

    @ConfigOption
    @Description("""
        One or more commands to be added to the Dockerfile used to build a Conda based image.
    """)
    final List<String> commands

    WaveBuildCondaConfig(Map opts) {
        mambaImage = opts.mambaImage
        baseImage = opts.baseImage
        basePackages = opts.basePackages
        commands = opts.commands as List<String>
    }

    /**
     * Adapt this config scope to the wave-api {@link CondaOpts} domain object,
     * applying the same defaults as {@link CondaOpts} for options that are not set.
     */
    CondaOpts toCondaOpts() {
        final Map<String,Object> opts = [:]
        if( mambaImage != null )
            opts.mambaImage = mambaImage
        if( baseImage != null )
            opts.baseImage = baseImage
        if( basePackages != null )
            opts.basePackages = basePackages
        if( commands != null )
            opts.commands = commands
        return new CondaOpts(opts)
    }
}


/**
 * Model the `wave.build.compression` config scope.
 *
 * This is a Nextflow-side mirror of the wave-api {@link BuildCompression} type,
 * which cannot implement {@link ConfigScope} because it belongs to an external
 * library. Use {@link #toBuildCompression()} to obtain the corresponding domain
 * object.
 */
@ToString(includeNames = true, includePackage = false, includeFields = true)
@CompileStatic
class WaveBuildCompressionConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        The compression algorithm used when building the container. Allowed values are `gzip`, `estargz` and `zstd` (default: `gzip`).
    """)
    final String mode

    @ConfigOption
    @Description("""
        The compression level used when building a container, depending on the chosen algorithm: gzip and estargz (0-9), zstd (0-22).
    """)
    final Integer level

    @ConfigOption
    @Description("""
        Forcefully apply the compression option to all layers, including already existing layers (default: `false`).
    """)
    final Boolean force

    WaveBuildCompressionConfig(Map opts) {
        mode = opts?.mode
        level = opts?.level as Integer
        force = opts?.force as Boolean
    }

    /**
     * Adapt this config scope to the wave-api {@link BuildCompression} domain
     * object, or {@code null} if no compression options are set.
     */
    BuildCompression toBuildCompression() {
        if( !mode && !level && !force )
            return null
        final result = new BuildCompression()
        if( mode )
            result.mode = BuildCompression.Mode.valueOf(mode.toLowerCase())
        if( level )
            result.level = level
        if( force )
            result.force = force
        return result
    }
}
