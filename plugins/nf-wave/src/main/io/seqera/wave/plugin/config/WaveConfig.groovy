/*
 * Copyright 2013-2023, Seqera Labs
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
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import io.seqera.wave.config.CondaOpts
import io.seqera.wave.config.SpackOpts
import nextflow.file.FileHelper
import nextflow.util.Duration
/**
 * Model Wave client configuration
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includeNames = true, includePackage = false, includeFields = true, excludes = 'reportOpts')
@CompileStatic
class WaveConfig {
    final private static String DEF_ENDPOINT = 'https://wave.seqera.io'
    final private static List<String> DEF_STRATEGIES = List.of('container','dockerfile','conda', 'spack')
    final private Boolean enabled
    final private String endpoint
    final private List<URL> containerConfigUrl
    final private Duration tokensCacheMaxDuration
    final private CondaOpts condaOpts
    final private SpackOpts spackOpts
    final private List<String> strategy
    final private Boolean bundleProjectResources
    final private String buildRepository
    final private String cacheRepository
    final private ReportOpts reportOpts
    final private RetryOpts retryOpts
    final private HttpOpts httpClientOpts
    final private Boolean freezeMode

    WaveConfig(Map opts, Map<String,String> env=System.getenv()) {
        this.enabled = opts.enabled
        this.endpoint = (opts.endpoint?.toString() ?: env.get('WAVE_API_ENDPOINT') ?: DEF_ENDPOINT)?.stripEnd('/')
        this.freezeMode = opts.freeze as Boolean
        this.containerConfigUrl = parseConfig(opts, env)
        this.tokensCacheMaxDuration = opts.navigate('tokens.cache.maxDuration', '30m') as Duration
        this.condaOpts = opts.navigate('build.conda', Collections.emptyMap()) as CondaOpts
        this.spackOpts = opts.navigate('build.spack', Collections.emptyMap()) as SpackOpts
        this.buildRepository = opts.navigate('build.repository') as String
        this.cacheRepository = opts.navigate('build.cacheRepository') as String
        this.strategy = parseStrategy(opts.strategy)
        this.bundleProjectResources = opts.bundleProjectResources
        this.reportOpts = new ReportOpts(opts.report as Map ?: Map.of())
        this.retryOpts = retryOpts0(opts)
        this.httpClientOpts = new HttpOpts(opts.httpClient as Map ?: Map.of())
        // some validation
        validateConfig()
    }

    Boolean enabled() { this.enabled }

    String endpoint() { this.endpoint }

    CondaOpts condaOpts() { this.condaOpts }

    SpackOpts spackOpts() { this.spackOpts }

    RetryOpts retryOpts() { this.retryOpts }

    HttpOpts httpOpts() { this.httpClientOpts }

    List<String> strategy() { this.strategy }

    boolean freezeMode() { return this.freezeMode }

    boolean bundleProjectResources() { bundleProjectResources }

    String buildRepository() { buildRepository }

    String cacheRepository() { cacheRepository }

    private void validateConfig() {
        def scheme= FileHelper.getUrlProtocol(endpoint)
        if( scheme !in ['http','https'] )
            throw new IllegalArgumentException("Endpoint URL should start with 'http:' or 'https:' protocol prefix - offending value: '$endpoint'")
        if( FileHelper.getUrlProtocol(buildRepository) )
            throw new IllegalArgumentException("Config setting 'wave.build.repository' should not include any protocol prefix - offending value: '$buildRepository'")
        if( FileHelper.getUrlProtocol(cacheRepository) )
            throw new IllegalArgumentException("Config setting 'wave.build.cacheRepository' should not include any protocol prefix - offending value: '$cacheRepository'")
    }

    private RetryOpts retryOpts0(Map opts) {
        if( opts.retryPolicy )
            return new RetryOpts(opts.retryPolicy as Map)
        if( opts.retry ) {
            log.warn "Configuration options 'wave.retry' has been deprecated - replace it with 'wave.retryPolicy'"
            return new RetryOpts(opts.retry as Map)
        }
        return new RetryOpts(Map.of())
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

    @Deprecated
    ReportOpts reportOpts() { reportOpts }
}
