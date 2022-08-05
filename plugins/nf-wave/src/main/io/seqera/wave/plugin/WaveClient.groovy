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

import javax.annotation.Nullable
import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.BasicFileAttributes
import java.time.Duration
import java.util.concurrent.Callable
import java.util.concurrent.TimeUnit

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.seqera.wave.plugin.util.DigestFunctions
import nextflow.Session
import nextflow.script.bundle.ModuleBundle
import nextflow.util.CacheHelper
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Wave client service
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WaveClient {

    private static Logger log = LoggerFactory.getLogger(WaveClient)

    final private HttpClient httpClient

    final private WaveConfig config

    final private String endpoint

    private Cache<String, SubmitContainerTokenResponse> cache

    WaveClient(Session session) {
        this.config = new WaveConfig(session.config.wave as Map ?: [:])
        this.endpoint = config.endpoint()
        log.debug "Wave server endpoint: ${endpoint}"
        // create cache
        cache = CacheBuilder<String, SubmitContainerTokenResponse>
            .newBuilder()
            .expireAfterWrite(config.tokensCacheMaxDuration().toSeconds(), TimeUnit.SECONDS)
            .build()
        // create http client
        this.httpClient = HttpClient.newBuilder()
                .version(HttpClient.Version.HTTP_1_1)
                .followRedirects(HttpClient.Redirect.NORMAL)
                .connectTimeout(Duration.ofSeconds(10))
                .build()
    }

    WaveConfig config() { return config }

    Boolean enabled() { config.enabled() }

    protected <T extends OutputStream> T makeTar(ModuleBundle bundle, T target) {
        try ( final archive = new TarArchiveOutputStream(target) ) {

            for (String name : bundle.getEntries() ) {
                final targetPath = bundle.path(name)
                final attrs = Files.readAttributes(targetPath, BasicFileAttributes)
                final entry = new TarArchiveEntry(targetPath, name)
                entry.setIds(0,0)
                entry.setGroupName("root")
                entry.setUserName("root")
                entry.setModTime(attrs.lastModifiedTime())
                entry.setMode(getMode(targetPath))
                // file permissions
                archive.putArchiveEntry(entry)
                if( !targetPath.isDirectory()) {
                    Files.copy(targetPath, archive)
                }
                archive.closeArchiveEntry()
            }
            archive.finish()
        }

        return target
    }

    /**
     * See {@link TarArchiveEntry#DEFAULT_DIR_MODE}
     */
    private static final int DIR_MODE = 040000;

    /**
     * See {@link TarArchiveEntry#DEFAULT_FILE_MODE}
     */
    private static final int FILE_MODE = 0100000;

    private int getMode(Path path) {
        final mode = path.isDirectory() ? DIR_MODE : FILE_MODE
        return mode + path.getPermissionsMode()
    }

    protected <T extends OutputStream> T makeGzip(InputStream source, T target) {
        try (final compressed = new GzipCompressorOutputStream(target)) {
            source.transferTo(compressed)
            compressed.flush()
        }
        return target
    }

    protected ContainerLayer layer(ModuleBundle bundle) {
        final tar = makeTar(bundle, new ByteArrayOutputStream()).toByteArray()
        final tarDigest = DigestFunctions.digest(tar)
        final gzip = makeGzip(new ByteArrayInputStream(tar), new ByteArrayOutputStream()).toByteArray()
        final gzipSize = gzip.length
        final gzipDigest = DigestFunctions.digest(gzip)
        final data = 'data:' + gzip.encodeBase64()

        log.debug """\
            Module bundle: ${bundle.root}
            - digest     : ${bundle.fingerprint()}    
            - location   : $data
            - tar digest : $tarDigest
            - gzip digest: $gzipDigest
            - gzip size  : $gzipSize
            """.stripIndent().rightTrim()

        return new ContainerLayer(
                location: data,
                tarDigest: tarDigest,
                gzipSize: gzipSize,
                gzipDigest: gzipDigest )

    }

    SubmitContainerTokenRequest makeRequest(ModuleBundle bundle, @Nullable String container, @Nullable ContainerConfig config) {
        if( config == null ) {
            config = new ContainerConfig()
        }
        // pre-prepend the bundle layer
        if( bundle && bundle.hasEntries() ) {
            config.prependLayer(layer(bundle))
        }
        final dockerContent = bundle?.dockerfile?.text?.bytes?.encodeBase64()?.toString()
        return new SubmitContainerTokenRequest(containerImage: container, containerConfig: config, containerFile: dockerContent)
    }

    SubmitContainerTokenResponse sendRequest(ModuleBundle bundle, @Nullable String container, @Nullable ContainerConfig config) {
        final req = makeRequest(bundle, container, config)
        return sendRequest(req)
    }

    SubmitContainerTokenResponse sendRequest(String image) {
        final configUrl = config().containerConfigUrl()
        final ContainerConfig containerConfig = configUrl ? fetchContainerConfig(configUrl) : null
        final request = new SubmitContainerTokenRequest(containerImage: image, containerConfig: containerConfig)
        return sendRequest(request)
    }

    SubmitContainerTokenResponse sendRequest(SubmitContainerTokenRequest request) {
        assert endpoint, 'Missing wave endpoint'
        assert !endpoint.endsWith('/'), "Endpoint url must not end with a slash - offending value: $endpoint"

        final body = JsonOutput.toJson(request)
        final uri = URI.create("${endpoint}/container-token")
        log.debug "Wave request: $uri - request: $request"
        final req = HttpRequest.newBuilder()
                .uri(uri)
                .headers('Content-Type','application/json')
                .POST(HttpRequest.BodyPublishers.ofString(body))
                .build()

        final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
        if( resp.statusCode()==200 ) {
            log.debug "Wave response: ${resp.body()}"
            return new JsonSlurper().parseText(resp.body()) as SubmitContainerTokenResponse
        }
        throw new BadResponseException("Wave invalid response: [${resp.statusCode()}] ${resp.body()}")
    }

    @Memoized
    ContainerConfig fetchContainerConfig(URL configUrl) {
        log.debug "Wave request container config: $configUrl"
        final req = HttpRequest.newBuilder()
                .uri(configUrl.toURI())
                .headers('Content-Type','application/json')
                .GET()
                .build()

        final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
        if( resp.statusCode()==200 ) {
            log.debug "Wave container config response: ${resp.body()}"
            return new JsonSlurper().parseText(resp.body()) as ContainerConfig
        }
        else {
            log.warn "Wave container config error response: [${resp.statusCode()}] ${resp.body()}"
            return null
        }
    }

    @Memoized
    protected String key0(ModuleBundle bundle, ContainerConfig containerConfig, String image) {
        final allMeta = new ArrayList(10)
        allMeta.add( bundle?.fingerprint() )
        allMeta.add( containerConfig?.hashCode() )
        allMeta.add( image )
        return CacheHelper.hasher(allMeta).hash().toString()
    }

    String fetchContainerImage(ModuleBundle bundle, String container, URL configUrl) {
        final ContainerConfig containerConfig = configUrl ? fetchContainerConfig(configUrl) : null
        // go ahead
        final key = key0(bundle, containerConfig, container)
        final result = cache.get(key, { sendRequest(bundle, container, containerConfig) } as Callable )
        return result.targetImage
    }

}
