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
import java.time.OffsetDateTime
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
import io.seqera.wave.plugin.config.FusionConfig
import io.seqera.wave.plugin.config.TowerConfig
import io.seqera.wave.plugin.config.WaveConfig
import io.seqera.wave.plugin.exception.BadResponseException
import io.seqera.wave.plugin.exception.UnauthorizedException
import io.seqera.wave.plugin.packer.Packer
import nextflow.Session
import nextflow.SysEnv
import nextflow.container.resolver.ContainerInfo
import nextflow.processor.TaskRun
import nextflow.script.bundle.ResourcesBundle
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

    private static final List<String> DEFAULT_CONDA_CHANNELS = ['conda-forge','defaults']

    final private HttpClient httpClient

    final private WaveConfig config

    final private FusionConfig fusion

    final private TowerConfig tower

    final private Packer packer

    final private String endpoint

    private Cache<String, SubmitContainerTokenResponse> cache

    private Session session

    private volatile String accessToken

    private volatile String refreshToken

    private CookieManager cookieManager

    private List<String> condaChannels

    WaveClient(Session session) {
        this.session = session
        this.config = new WaveConfig(session.config.wave as Map ?: Collections.emptyMap(), SysEnv.get())
        this.fusion = new FusionConfig(session.config.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
        this.tower = new TowerConfig(session.config.tower as Map ?: Collections.emptyMap(), SysEnv.get())
        this.endpoint = config.endpoint()
        this.condaChannels = session.getCondaConfig()?.getChannels() ?: DEFAULT_CONDA_CHANNELS
        log.debug "Wave server endpoint: ${endpoint}"
        this.packer = new Packer()
        // create cache
        cache = CacheBuilder<String, SubmitContainerTokenResponse>
            .newBuilder()
            .expireAfterWrite(config.tokensCacheMaxDuration().toSeconds(), TimeUnit.SECONDS)
            .build()
        // the cookie manager
        cookieManager = new CookieManager()
        // create http client
        this.httpClient = HttpClient.newBuilder()
                .version(HttpClient.Version.HTTP_1_1)
                .followRedirects(HttpClient.Redirect.NORMAL)
                .cookieHandler(cookieManager)
                .connectTimeout(Duration.ofSeconds(10))
                .build()
    }

    WaveConfig config() { return config }

    Boolean enabled() { return config.enabled() }

    protected ContainerLayer makeLayer(ResourcesBundle bundle) {
        final result = packer.layer(bundle.content())
        return result
    }

    SubmitContainerTokenRequest makeRequest(WaveAssets assets) {
        final containerConfig = assets.containerConfig ?: new ContainerConfig()
        // prepend the bundle layer
        if( assets.moduleResources!=null && assets.moduleResources.hasEntries() ) {
            containerConfig.prependLayer(makeLayer(assets.moduleResources))
        }
        // prepend project resources bundle
        if( assets.projectResources!=null && assets.projectResources.hasEntries() ) {
            containerConfig.prependLayer(makeLayer(assets.projectResources))
        }

        if( !assets.containerImage && !assets.dockerFileContent )
            throw new IllegalArgumentException("Wave container request requires at least a image or container file to build")

        if( assets.containerImage && assets.dockerFileContent )
            throw new IllegalArgumentException("Wave container image and container file cannot be specified in the same request")

        return new SubmitContainerTokenRequest(
                containerImage: assets.containerImage,
                containerPlatform: assets.containerPlatform,
                containerConfig: containerConfig,
                containerFile: assets.dockerFileEncoded(),
                condaFile: assets.condaFileEncoded(),
                buildRepository: config().buildRepository(),
                cacheRepository: config.cacheRepository(),
                timestamp: OffsetDateTime.now().toString(),
                fingerprint: assets.fingerprint()
        )
    }

    SubmitContainerTokenResponse sendRequest(WaveAssets assets) {
        final req = makeRequest(assets)
        req.towerAccessToken = tower.accessToken
        req.towerWorkspaceId = tower.workspaceId
        req.towerEndpoint = tower.endpoint
        return sendRequest(req)
    }

    SubmitContainerTokenResponse sendRequest(String image) {
        final ContainerConfig containerConfig = resolveContainerConfig()
        final request = new SubmitContainerTokenRequest(
                containerImage: image,
                containerConfig: containerConfig,
                towerAccessToken: tower.accessToken,
                towerWorkspaceId: tower.workspaceId,
                towerEndpoint: tower.endpoint )
        return sendRequest(request)
    }

    SubmitContainerTokenResponse sendRequest(SubmitContainerTokenRequest request) {
        return sendRequest0(request, 1)
    }

    SubmitContainerTokenResponse sendRequest0(SubmitContainerTokenRequest request, int attempt) {
        assert endpoint, 'Missing wave endpoint'
        assert !endpoint.endsWith('/'), "Endpoint url must not end with a slash - offending value: $endpoint"

        // update the request token
        accessToken ?= tower.accessToken
        refreshToken ?= tower.refreshToken

        // set the request access token
        request.towerAccessToken = accessToken

        final body = JsonOutput.toJson(request)
        final uri = URI.create("${endpoint}/container-token")
        log.debug "Wave request: $uri; attempt=$attempt - request: $request"
        final req = HttpRequest.newBuilder()
                .uri(uri)
                .headers('Content-Type','application/json')
                .POST(HttpRequest.BodyPublishers.ofString(body))
                .build()

        try {
            final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
            log.debug "Wave response: statusCode=${resp.statusCode()}; body=${resp.body()}"
            if( resp.statusCode()==200 )
                return jsonToSubmitResponse(resp.body())
            if( resp.statusCode()==401 ) {
                final shouldRetry = request.towerAccessToken
                        && refreshToken
                        && attempt==1
                        && refreshJwtToken0(refreshToken)
                if( shouldRetry ) {
                    return sendRequest0(request, attempt+1)
                }
                else
                    throw new UnauthorizedException("Unauthorised [401] - Verify you have provided a valid access token")
            }
            else
                throw new BadResponseException("Wave invalid response: [${resp.statusCode()}] ${resp.body()}")
        }
        catch (ConnectException e) {
            throw new IllegalStateException("Unable to connect Wave service: $endpoint")
        }
    }

    private SubmitContainerTokenResponse jsonToSubmitResponse(String body) {
        final type = new TypeToken<SubmitContainerTokenResponse>(){}.getType()
        return new Gson().fromJson(body, type)
    }

    private ContainerConfig jsonToContainerConfig(String json) {
        final type = new TypeToken<ContainerConfig>(){}.getType()
        return new Gson().fromJson(json, type)
    }

    protected URL defaultFusionUrl() {
        final isArm = config.containerPlatform()?.tokenize('/')?.contains('arm64')
        return isArm
                ? new URL(FusionConfig.DEFAULT_FUSION_ARM64_URL)
                : new URL(FusionConfig.DEFAULT_FUSION_AMD64_URL)
    }

    ContainerConfig resolveContainerConfig() {
        final urls = new ArrayList<URL>(config.containerConfigUrl())
        if( fusion.enabled() ) {
            final fusionUrl = fusion.containerConfigUrl() ?: defaultFusionUrl()
            urls.add(fusionUrl)
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
    synchronized protected ContainerConfig fetchContainerConfig(URL configUrl) {
        log.debug "Wave request container config: $configUrl"
        final req = HttpRequest.newBuilder()
                .uri(configUrl.toURI())
                .headers('Content-Type','application/json')
                .GET()
                .build()

        final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
        final code = resp.statusCode()
        if( code>=200 && code<400 ) {
            log.debug "Wave container config response: [$code] ${resp.body()}"
            return jsonToContainerConfig(resp.body())
        }
        throw new BadResponseException("Unexpected response for containerContainerConfigUrl \'$configUrl\': [${resp.statusCode()}] ${resp.body()}")
    }

    protected void checkConflicts(Map<String,String> attrs, String name) {
        if( attrs.dockerfile && attrs.conda ) {
            throw new IllegalArgumentException("Process '${name}' declares both a 'conda' directive and a module bundle dockerfile that conflicts each other")
        }
        if( attrs.container && attrs.dockerfile ) {
            throw new IllegalArgumentException("Process '${name}' declares both a 'container' directive and a module bundle dockerfile that conflicts each other")
        }
        if( attrs.container && attrs.conda ) {
            throw new IllegalArgumentException("Process '${name}' declares both 'container' and 'conda' directives that conflicts each other")
        }
    }

    Map<String,String> resolveConflicts(Map<String,String> attrs, List<String> strategy) {
        final result = new HashMap<String,String>()
        for( String it : strategy ) {
            if( attrs.get(it) ) {
                return [(it): attrs.get(it)]
            }
        }
        return result
    }

    @Memoized
    WaveAssets resolveAssets(TaskRun task, String containerImage) {
        // get the bundle
        final bundle = task.getModuleBundle()
        // compose the request attributes
        def attrs = new HashMap<String,String>()
        attrs.container = containerImage
        attrs.conda = task.config.conda as String
        if( bundle!=null && bundle.dockerfile ) {
            attrs.dockerfile = bundle.dockerfile.text
        }

        // validate request attributes
        if( config().strategy() )
            attrs = resolveConflicts(attrs, config().strategy())
        else
            checkConflicts(attrs, task.lazyName())

        //  resolve the wave assets
        return resolveAssets0(attrs, bundle)
    }

    protected WaveAssets resolveAssets0(Map<String,String> attrs, ResourcesBundle bundle) {

        String dockerScript = attrs.dockerfile
        final containerImage = attrs.container

        /*
         * If 'conda' directive is specified use it to create a Dockefile
         * to assemble the target container
         */
        Path condaFile = null
        if( attrs.conda ) {
            if( dockerScript )
                throw new IllegalArgumentException("Unexpected conda and dockerfile conflict")

            // map the recipe to a dockerfile
            if( isCondaFile(attrs.conda) ) {
                condaFile = Path.of(attrs.conda)
                dockerScript = condaFileToDockerFile()
            }
            else {
                dockerScript = condaRecipeToDockerFile(attrs.conda)
            }
        }

        /*
         * The process should declare at least a container image name via 'container' directive
         * or a dockerfile file to build, otherwise there's no job to be done by wave
         */
        if( !dockerScript && !containerImage ) {
            return null
        }

        /*
         * project level resources i.e. `ROOT/bin/` directory files
         * are only uploaded when using fusion
         */
        final projectRes = config.bundleProjectResources() && session.binDir
                    ? projectResources(session.binDir)
                    : null

        /*
         * the container platform to be used
         */
        final platform = config.containerPlatform()

        // read the container config and go ahead
        final containerConfig = this.resolveContainerConfig()
        return new WaveAssets(
                    containerImage,
                    platform,
                    bundle,
                    containerConfig,
                    dockerScript,
                    condaFile,
                    projectRes)
    }

    @Memoized
    protected ResourcesBundle projectResources(Path path) {
        log.debug "Wave assets bundle bin resources: $path"
        // place project 'bin' resources under '/usr/local'
        // see https://unix.stackexchange.com/questions/8656/usr-bin-vs-usr-local-bin-on-linux
        return path && path.parent
                ? ResourcesBundle.scan( path.parent, [filePattern: "$path.name/**", baseDirectory:'usr/local'] )
                : null
    }

    ContainerInfo fetchContainerImage(WaveAssets assets) {
        try {
            // compute a unique hash for this request assets
            final key = assets.fingerprint()
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
        FROM ${config.condaOpts().mambaImage}
        COPY --chown=\$MAMBA_USER:\$MAMBA_USER conda.yml /tmp/conda.yml
        RUN micromamba install -y -n base -f /tmp/conda.yml && \\
            micromamba clean -a -y
        """.stripIndent()

        return addCommands(result)
    }

    protected String addCommands(String result) {
        if( !config.condaOpts().commands )
            return result
        for( String cmd : config.condaOpts().commands ) {
            result += cmd + "\n"
        }
        return result
    }

    protected String condaRecipeToDockerFile(String recipe) {
        def channelsOpts = condaChannels.collect(it -> "-c $it").join(' ')
        def result = """\
        FROM ${config.condaOpts().mambaImage}
        RUN \\
           micromamba install -y -n base $channelsOpts \\
           $recipe \\
           && micromamba clean -a -y
        """.stripIndent()

        return addCommands(result)
    }

    protected boolean isCondaFile(String value) {
        if( value.contains('\n') )
            return false
        return value.endsWith('.yaml') || value.endsWith('.yml') || value.endsWith('.txt')
    }

    protected boolean refreshJwtToken0(String refresh) {
        log.debug "Token refresh request >> $refresh"

        final req = HttpRequest.newBuilder()
                .uri(new URI("${tower.endpoint}/oauth/access_token"))
                .headers('Content-Type',"application/x-www-form-urlencoded")
                .POST(HttpRequest.BodyPublishers.ofString("grant_type=refresh_token&refresh_token=${URLEncoder.encode(refresh, 'UTF-8')}"))
                .build()

        final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
        log.debug "Refresh cookie response: [${resp.statusCode()}] ${resp.body()}"
        if( resp.statusCode() != 200 )
            return false

        final authCookie = getCookie('JWT')
        final refreshCookie = getCookie('JWT_REFRESH_TOKEN')

        // set the new bearer token in the current client session
        if( authCookie?.value ) {
            log.trace "Updating http client bearer token=$authCookie.value"
            accessToken = authCookie.value
        }
        else {
            log.warn "Missing JWT cookie from refresh token response ~ $authCookie"
        }

        // set the new refresh token
        if( refreshCookie?.value ) {
            log.trace "Updating http client refresh token=$refreshCookie.value"
            refreshToken = refreshCookie.value
        }
        else {
            log.warn "Missing JWT_REFRESH_TOKEN cookie from refresh token response ~ $refreshCookie"
        }

        return true
    }

    private HttpCookie getCookie(final String cookieName) {
        for( HttpCookie it : cookieManager.cookieStore.cookies ) {
            if( it.name == cookieName )
                return it
        }
        return null
    }

}
