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

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Path
import java.time.Duration
import java.util.concurrent.Callable
import java.util.concurrent.TimeUnit

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.util.concurrent.UncheckedExecutionException
import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.seqera.wave.plugin.config.WaveConfig
import io.seqera.wave.plugin.packer.Packer
import nextflow.Session
import nextflow.container.resolver.ContainerInfo
import nextflow.processor.TaskRun
import nextflow.script.bundle.ModuleBundle
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

    final private FusionConfig fusion

    final private Packer packer

    final private String endpoint

    private Cache<String, SubmitContainerTokenResponse> cache

    WaveClient(Session session) {
        this.config = new WaveConfig(session.config.wave as Map ?: [:])
        this.fusion = new FusionConfig(session.config.fusion as Map ?: [:])
        this.endpoint = config.endpoint()
        log.debug "Wave server endpoint: ${endpoint}"
        this.packer = new Packer()
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
        packer.makeTar(bundle.root, bundle.getPathsList(), target)
    }

    protected ContainerLayer makeLayer(ModuleBundle bundle) {

        final result = packer.layer(bundle.root, bundle.pathsList)

        log.debug """\
            Module bundle: ${bundle.root}
            - digest     : ${bundle.fingerprint()}    
            - location   : ${result.location}
            - tar digest : ${result.tarDigest}
            - gzip digest: ${result.gzipDigest}
            - gzip size  : ${result.gzipSize}
            """.stripIndent().rightTrim()

        return result
    }

    SubmitContainerTokenRequest makeRequest(WaveAssets assets) {
        final containerConfig = assets.containerConfig ?: new ContainerConfig()
        // pre-prepend the bundle layer
        if( assets.bundle && assets.bundle.hasEntries() ) {
            containerConfig.prependLayer(makeLayer(assets.bundle))
        }

        if( !assets.containerImage && !assets.dockerFileContent )
            throw new IllegalArgumentException("Wave container request requires at least a image or container file to build")

        return new SubmitContainerTokenRequest(
                containerImage: assets.containerImage,
                containerConfig: containerConfig,
                containerFile: assets.dockerFileEncoded(),
                condaFile: assets.condaFileEncoded()
        )
    }

    SubmitContainerTokenResponse sendRequest(WaveAssets assets) {
        final req = makeRequest(assets)
        return sendRequest(req)
    }

    SubmitContainerTokenResponse sendRequest(String image) {
        final ContainerConfig containerConfig = resolveContainerConfig()
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
            return jsonToSubmitResponse(resp.body())
        }
        throw new BadResponseException("Wave invalid response: [${resp.statusCode()}] ${resp.body()}")
    }

    private SubmitContainerTokenResponse jsonToSubmitResponse(String body) {
        final type = new TypeToken<SubmitContainerTokenResponse>(){}.getType()
        return new Gson().fromJson(body, type)
    }

    private ContainerConfig jsonToContainerConfig(String json) {
        final type = new TypeToken<ContainerConfig>(){}.getType()
        return new Gson().fromJson(json, type)
    }

    ContainerConfig resolveContainerConfig() {
        final urls = new ArrayList<URL>(config.containerConfigUrl())
        if( fusion.enabled() && fusion.containerConfigUrl() ) {
            urls.add( fusion.containerConfigUrl() )
        }
        if( !urls )
            return null
        def result = new ContainerConfig()
        for( URL it : urls ) {
            // append each config to the other - the last has priority
            result += fetchContainerConfig(it)
        }
        return result
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
            return jsonToContainerConfig(resp.body())
        }
        else {
            log.warn "Wave container config error response: [${resp.statusCode()}] ${resp.body()}"
            return null
        }
    }

    @Memoized
    WaveAssets resolveAssets(TaskRun task, String containerImage) {
        final bundle = task.getModuleBundle()
        // read the dockerfile defined in the module bundle
        String dockerFileContent = null
        if( bundle!=null && bundle.dockerfile ) {
            dockerFileContent = bundle.dockerfile.text
        }
        // compute docker file + conda content
        final condaRecipe = task.config.conda as String
        if( dockerFileContent && condaRecipe ) {
            throw new IllegalArgumentException("Process '${task.lazyName()}' declares both a 'conda' directive and a module bundle dockerfile that conflicts each other")
        }
        Path condaFile = null
        if( condaRecipe ) {
            // map the recipe to a dockerfile
            if( isCondaFile(condaRecipe) ) {
                condaFile = Path.of(condaRecipe)
                dockerFileContent = condaFileToDockerFile()
            }
            else {
                dockerFileContent = condaRecipeToDockerFile(condaRecipe)
            }
        }

        if( dockerFileContent && containerImage ) {
            throw new IllegalArgumentException("Process '${task.lazyName()}' declares both a 'container' directive and a module bundle dockerfile that conflicts each other")
        }

        if( !dockerFileContent && !containerImage ) {
            // nothing to do
            return null
        }

        // read the container config and go ahead
        final containerConfig = this.resolveContainerConfig()
        return new WaveAssets(containerImage, bundle, containerConfig, dockerFileContent, condaFile)
    }

    ContainerInfo fetchContainerImage(WaveAssets assets) {
        try {
            // compute a unique hash for this request assets
            final key = assets.hashKey()
            // get from cache or submit a new request
            final response = cache.get(key, { sendRequest(assets) } as Callable )
            // assemble the container info response
            return new ContainerInfo(assets.containerImage, response.targetImage, key)
        }
        catch ( UncheckedExecutionException e ) {
            throw e.cause
        }
    }

    protected String condaFileToDockerFile() {
        def result = """\
        FROM ${config.mambaOpts().from}
        COPY --chown=\$MAMBA_USER:\$MAMBA_USER conda.yml /tmp/conda.yml
        RUN micromamba install -y -n base -f /tmp/conda.yml && \\
            micromamba clean -a -y
        """.stripIndent()

        return addUser(result)
    }

    protected String addUser(String result) {
        if( config.mambaOpts().user )
            result += "USER ${config.mambaOpts().user}\n"
        return result
    }

    protected String condaRecipeToDockerFile(String recipe) {
        def result = """\
        FROM ${config.mambaOpts().from}
        RUN \\
           micromamba install -y -n base -c defaults -c conda-forge \\
           $recipe \\
           && micromamba clean -a -y
        """.stripIndent()

        return addUser(result)
    }

    protected boolean isCondaFile(String value) {
        if( value.contains('\n') )
            return false
        return value.endsWith('.yaml') || value.endsWith('.yml') || value.endsWith('.txt')
    }
}
