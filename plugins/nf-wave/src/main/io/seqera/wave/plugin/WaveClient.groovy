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

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import nextflow.Session
import io.seqera.wave.plugin.util.DigestFunctions
import nextflow.script.bundle.ModuleBundle
import org.apache.commons.compress.archivers.tar.TarArchiveEntry
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class WaveClient {

    private static Logger log = LoggerFactory.getLogger(WaveClient)

    final private HttpClient httpClient

    final private WaveConfig config

    final private String endpoint

    WaveClient(Session session) {
        this.config = new WaveConfig(session.config.wave as Map)
        this.endpoint = config.endpoint()
        log.debug "Wave server endpoint: ${endpoint}"
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

    SubmitContainerTokenRequest makeRequest(ModuleBundle bundle, @Nullable String container) {
        final layers = bundle ? [layer(bundle)] : []
        final config = new ContainerConfig(layers: layers)
        final dockerContent = bundle.dockerfile?.text?.bytes?.encodeBase64()?.toString()
        return new SubmitContainerTokenRequest(containerImage: container, containerConfig: config, containerFile: dockerContent)
    }

    SubmitContainerTokenResponse sendRequest(ModuleBundle bundle, @Nullable String container) {
        final req = makeRequest(bundle, container)
        return sendRequest(req)
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
        else {
            log.warn "Wave error response: [${resp.statusCode()}] ${resp.body()}"
            return null
        }
    }
}
