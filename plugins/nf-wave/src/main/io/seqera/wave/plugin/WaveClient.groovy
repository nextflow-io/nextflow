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

package io.seqera.wave.plugin

import static io.seqera.wave.util.DockerHelper.*

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Path
import java.time.Duration
import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.function.Predicate

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.util.concurrent.UncheckedExecutionException
import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.seqera.wave.plugin.config.TowerConfig
import io.seqera.wave.plugin.config.WaveConfig
import io.seqera.wave.plugin.exception.BadResponseException
import io.seqera.wave.plugin.exception.UnauthorizedException
import io.seqera.wave.plugin.packer.Packer
import nextflow.Session
import nextflow.SysEnv
import nextflow.container.inspect.ContainerInspectMode
import nextflow.container.resolver.ContainerInfo
import nextflow.fusion.FusionConfig
import nextflow.processor.Architecture
import nextflow.processor.TaskRun
import nextflow.script.bundle.ResourcesBundle
import nextflow.util.SysHelper
import nextflow.util.Threads
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

    final static private String[] REQUEST_HEADERS =  new String[]{
                        'Content-Type','application/json',
                        'Accept','application/json',
                        'Accept','application/vnd.oci.image.index.v1+json',
                        'Accept','application/vnd.oci.image.manifest.v1+json',
                        'Accept','application/vnd.docker.distribution.manifest.v1+prettyjws',
                        'Accept','application/vnd.docker.distribution.manifest.v2+json',
                        'Accept','application/vnd.docker.distribution.manifest.list.v2+json' }

    private static final List<String> DEFAULT_CONDA_CHANNELS = ['seqera','conda-forge','bioconda','defaults']

    private static final String DEFAULT_SPACK_ARCH = 'x86_64'

    private static final String DEFAULT_DOCKER_PLATFORM = 'linux/amd64'

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

    final private String waveRegistry

    WaveClient(Session session) {
        this.session = session
        this.config = new WaveConfig(session.config.wave as Map ?: Collections.emptyMap(), SysEnv.get())
        this.fusion = new FusionConfig(session.config.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
        this.tower = new TowerConfig(session.config.tower as Map ?: Collections.emptyMap(), SysEnv.get())
        this.endpoint = config.endpoint()
        this.condaChannels = session.getCondaConfig()?.getChannels() ?: DEFAULT_CONDA_CHANNELS
        log.debug "Wave config: $config"
        this.packer = new Packer()
        this.waveRegistry = new URI(endpoint).getAuthority()
        // create cache
        cache = CacheBuilder<String, SubmitContainerTokenResponse>
            .newBuilder()
            .expireAfterWrite(config.tokensCacheMaxDuration().toSeconds(), TimeUnit.SECONDS)
            .build()
        // the cookie manager
        cookieManager = new CookieManager()
        // create http client
        this.httpClient = newHttpClient()
    }

    protected HttpClient newHttpClient() {
        final builder = HttpClient.newBuilder()
                .version(HttpClient.Version.HTTP_1_1)
                .followRedirects(HttpClient.Redirect.NEVER)
                .cookieHandler(cookieManager)
                .connectTimeout(config.httpOpts().connectTimeout())
        // use virtual threads executor if enabled
        if( Threads.useVirtual() )
            builder.executor(Executors.newVirtualThreadPerTaskExecutor())
        // build and return the new client
        return builder.build()
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

        if( !assets.containerImage && !assets.containerFile )
            throw new IllegalArgumentException("Wave container request requires at least a image or container file to build")

        if( assets.containerImage && assets.containerFile )
            throw new IllegalArgumentException("Wave container image and container file cannot be specified in the same request")

        return new SubmitContainerTokenRequest(
                containerImage: assets.containerImage,
                containerPlatform: assets.containerPlatform,
                containerConfig: containerConfig,
                containerFile: assets.dockerFileEncoded(),
                condaFile: assets.condaFileEncoded(),
                spackFile: assets.spackFileEncoded(),
                buildRepository: config().buildRepository(),
                cacheRepository: config.cacheRepository(),
                timestamp: OffsetDateTime.now().toString(),
                fingerprint: assets.fingerprint(),
                freeze: config.freezeMode(),
                format: assets.singularity ? 'sif' : null,
                dryRun: ContainerInspectMode.active()
        )
    }

    SubmitContainerTokenResponse sendRequest(WaveAssets assets) {
        final req = makeRequest(assets)
        req.towerAccessToken = tower.accessToken
        req.towerRefreshToken = tower.refreshToken
        req.towerWorkspaceId = tower.workspaceId
        req.towerEndpoint = tower.endpoint
        req.workflowId = tower.workflowId
        return sendRequest(req)
    }

    SubmitContainerTokenResponse sendRequest(String image) {
        final ContainerConfig containerConfig = resolveContainerConfig()
        final request = new SubmitContainerTokenRequest(
                containerImage: image,
                containerConfig: containerConfig,
                towerAccessToken: tower.accessToken,
                towerWorkspaceId: tower.workspaceId,
                towerEndpoint: tower.endpoint,
                workflowId: tower.workflowId,
                freeze: config.freezeMode(),
                dryRun: ContainerInspectMode.active(),
        )
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
        request.towerRefreshToken = refreshToken

        final body = JsonOutput.toJson(request)
        final uri = URI.create("${endpoint}/container-token")
        log.debug "Wave request: $uri; attempt=$attempt - request: $request"
        final req = HttpRequest.newBuilder()
                .uri(uri)
                .headers('Content-Type','application/json')
                .POST(HttpRequest.BodyPublishers.ofString(body))
                .build()

        try {
            final resp = httpSend(req)
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
                    throw new UnauthorizedException("Unauthorized [401] - Verify you have provided a valid access token")
            }
            else
                throw new BadResponseException("Wave invalid response: [${resp.statusCode()}] ${resp.body()}")
        }
        catch (IOException e) {
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

    protected URL defaultFusionUrl(String platform) {
        final isArm = platform.tokenize('/')?.contains('arm64')
        return isArm
                ? new URL(FusionConfig.DEFAULT_FUSION_ARM64_URL)
                : new URL(FusionConfig.DEFAULT_FUSION_AMD64_URL)
    }

    ContainerConfig resolveContainerConfig(String platform = DEFAULT_DOCKER_PLATFORM) {
        final urls = new ArrayList<URL>(config.containerConfigUrl())
        if( fusion.enabled() ) {
            final fusionUrl = fusion.containerConfigUrl() ?: defaultFusionUrl(platform)
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

        final resp = httpSend(req)
        final code = resp.statusCode()
        final body = resp.body()
        if( code>=200 && code<400 ) {
            log.debug "Wave container config response: [$code] ${body}"
            return jsonToContainerConfig(resp.body())
        }
        throw new BadResponseException("Unexpected response for containerContainerConfigUrl \'$configUrl\': [${resp.statusCode()}] ${resp.body()}")
    }

    protected void checkConflicts(Map<String,String> attrs, String name) {
        if( attrs.container && attrs.conda ) {
            throw new IllegalArgumentException("Process '${name}' declares both 'container' and 'conda' directives that conflict each other")
        }
        if( attrs.container && attrs.spack ) {
            throw new IllegalArgumentException("Process '${name}' declares both 'container' and 'spack' directives that conflict each other")
        }
        if( attrs.spack && attrs.conda ) {
            throw new IllegalArgumentException("Process '${name}' declares both 'spack' and 'conda' directives that conflict each other")
        }
        checkConflicts0(attrs, name, 'dockerfile')
        checkConflicts0(attrs, name, 'singularityfile')
    }

    protected void checkConflicts0(Map<String,String> attrs, String name, String fileType) {
        if( attrs.get(fileType) && attrs.conda ) {
            throw new IllegalArgumentException("Process '${name}' declares both a 'conda' directive and a module bundle $fileType that conflict each other")
        }
        if( attrs.container && attrs.get(fileType) ) {
            throw new IllegalArgumentException("Process '${name}' declares both a 'container' directive and a module bundle $fileType that conflict each other")
        }
        if( attrs.get(fileType) && attrs.spack ) {
            throw new IllegalArgumentException("Process '${name}' declares both a 'spack' directive and a module bundle $fileType that conflict each other")
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

    protected List<String> patchStrategy(List<String> strategy, boolean singularity) {
        if( !singularity )
            return strategy
        // when singularity is enabled, replaces `dockerfile` with `singularityfile`
        // in the strategy if not specified explicitly
        final p = strategy.indexOf('dockerfile')
        if( p!=-1 && !strategy.contains('singularityfile') ) {
            final result = new ArrayList(strategy)
            result.remove(p)
            result.add(p, 'singularityfile')
            return Collections.<String>unmodifiableList(result)
        }
        return strategy
    }

    static Architecture defaultArch() {
        try {
            return new Architecture(SysHelper.getArch())
        }
        catch (Exception e) {
            log.debug "Unable to detect system arch", e
            return null
        }
    }

    @Memoized
    WaveAssets resolveAssets(TaskRun task, String containerImage, boolean singularity) {
        // get the bundle
        final bundle = task.getModuleBundle()
        // get the Spack architecture
        final arch = task.config.getArchitecture()
        final spackArch = arch ? arch.spackArch : DEFAULT_SPACK_ARCH
        final dockerArch = arch? arch.dockerArch : DEFAULT_DOCKER_PLATFORM
        // compose the request attributes
        def attrs = new HashMap<String,String>()
        attrs.container = containerImage
        attrs.conda = task.config.conda as String
        attrs.spack = task.config.spack as String
        if( bundle!=null && bundle.dockerfile ) {
            attrs.dockerfile = bundle.dockerfile.text
        }
        if( bundle!=null && bundle.singularityfile ) {
            attrs.singularityfile = bundle.singularityfile.text
        }

        // validate request attributes
        final strategy = config().strategy()
        if( strategy )
            attrs = resolveConflicts(attrs, patchStrategy(strategy, singularity))
        else
            checkConflicts(attrs, task.lazyName())

        //  resolve the wave assets
        return resolveAssets0(attrs, bundle, singularity, dockerArch, spackArch)
    }

    protected WaveAssets resolveAssets0(Map<String,String> attrs, ResourcesBundle bundle, boolean singularity, String dockerArch, String spackArch) {

        final scriptType = singularity ? 'singularityfile' : 'dockerfile'
        String containerScript = attrs.get(scriptType)
        final containerImage = attrs.container

        /*
         * If 'conda' directive is specified use it to create a container file
         * to assemble the target container
         */
        Path condaFile = null
        if( attrs.conda ) {
            if( containerScript )
                throw new IllegalArgumentException("Unexpected conda and $scriptType conflict while resolving wave container")

            // map the recipe to a dockerfile
            if( isCondaRemoteFile(attrs.conda) ) {
                containerScript = singularity
                        ? condaPackagesToSingularityFile(attrs.conda, condaChannels, config.condaOpts())
                        : condaPackagesToDockerFile(attrs.conda, condaChannels, config.condaOpts())
            }
            else {
                if( isCondaLocalFile(attrs.conda) ) {
                    // 'conda' attribute is the path to the local conda environment
                    // note: ignore the 'channels' attribute because they are supposed to be provided by the conda file
                    condaFile = condaFileFromPath(attrs.conda, null)
                }
                else {
                    // 'conda' attributes is resolved as the conda packages to be used
                    condaFile = condaFileFromPackages(attrs.conda, condaChannels)
                }
                // create the container file to build the container
                containerScript = singularity
                        ? condaFileToSingularityFile(config.condaOpts())
                        : condaFileToDockerFile(config.condaOpts())
            }
        }

        /*
         * If 'spack' directive is specified use it to create a container file
         * to assemble the target container
         */
        Path spackFile = null
        if( attrs.spack ) {
            if( containerScript )
                throw new IllegalArgumentException("Unexpected spack and dockerfile conflict while resolving wave container")

            if( isSpackFile(attrs.spack) ) {
                // parse the attribute as a spack file path *and* append the base packages if any
                spackFile = addPackagesToSpackFile(attrs.spack, config.spackOpts())
            }
            else {
                // create a minimal spack file with package spec from user input
                spackFile = spackPackagesToSpackFile(attrs.spack, config.spackOpts())
            }
            // create the container file to build the container
            containerScript = singularity
                    ? spackFileToSingularityFile(config.spackOpts())
                    : spackFileToDockerFile(config.spackOpts())
        }

        /*
         * The process should declare at least a container image name via 'container' directive
         * or a dockerfile file to build, otherwise there's no job to be done by wave
         */
        if( !containerScript && !containerImage ) {
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
        final platform = dockerArch

        // check is a valid container image
        WaveAssets.validateContainerName(containerImage)
        // read the container config and go ahead
        final containerConfig = this.resolveContainerConfig(platform)
        return new WaveAssets(
                    containerImage,
                    platform,
                    bundle,
                    containerConfig,
                    containerScript,
                    condaFile,
                    spackFile,
                    projectRes,
                    singularity)
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
            if( config.freezeMode() )  {
                if( response.buildId && !ContainerInspectMode.active() ) {
                    // await the image to be available when a new image is being built
                    awaitImage(response.targetImage)
                }
                return new ContainerInfo(assets.containerImage, response.containerImage, key)
            }
            // assemble the container info response
            return new ContainerInfo(assets.containerImage, response.targetImage, key)
        }
        catch ( UncheckedExecutionException e ) {
            throw e.cause
        }
    }

    protected URI imageToManifestUri(String image) {
        final p = image.indexOf('/')
        if( p==-1 ) throw new IllegalArgumentException("Invalid container name: $image")
        final result = 'https://' + image.substring(0,p) + '/v2' + image.substring(p).replace(':','/manifests/')
        return new URI(result)
    }

    protected void awaitImage(String image) {
        final manifest = imageToManifestUri(image)
        final req = HttpRequest.newBuilder()
                .uri(manifest)
                .headers(REQUEST_HEADERS)
                .timeout(Duration.ofMinutes(5))
                .GET()
                .build()
        final begin = System.currentTimeMillis()
        final resp = httpSend(req)
        final body = resp.body()
        final code = resp.statusCode()
        if( code>=200 && code<400 ) {
            log.debug "Wave container available in ${nextflow.util.Duration.of(System.currentTimeMillis()-begin)}: [$code] ${body}"
        }
        else
            throw new BadResponseException("Unexpected response for \'$manifest\': [${code}] ${body}")
    }

    static protected boolean isCondaLocalFile(String value) {
        if( value.contains('\n') )
            return false
        if( value.startsWith('http://') || value.startsWith('https://') )
            return false
        return value.endsWith('.yaml') || value.endsWith('.yml') || value.endsWith('.txt')
    }

    static protected boolean isCondaRemoteFile(String value) {
        value.startsWith('http://') || value.startsWith('https://')
    }

    protected boolean isSpackFile(String value) {
        if( value.contains('\n') )
            return false
        return value.endsWith('.yaml') || value.endsWith('.yml')
    }

    protected boolean refreshJwtToken0(String refresh) {
        log.debug "Token refresh request >> $refresh"

        final req = HttpRequest.newBuilder()
                .uri(new URI("${tower.endpoint}/oauth/access_token"))
                .headers('Content-Type',"application/x-www-form-urlencoded")
                .POST(HttpRequest.BodyPublishers.ofString("grant_type=refresh_token&refresh_token=${URLEncoder.encode(refresh, 'UTF-8')}"))
                .build()

        final resp = httpSend(req)
        final code = resp.statusCode()
        final body = resp.body()
        log.debug "Refresh cookie response: [${code}] ${body}"
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

    protected <T> RetryPolicy<T> retryPolicy(Predicate<? extends Throwable> cond, Predicate<T> handle) {
        final cfg = config.retryOpts()
        final listener = new EventListener<ExecutionAttemptedEvent<T>>() {
            @Override
            void accept(ExecutionAttemptedEvent<T> event) throws Throwable {
                def msg = "Wave connection failure - attempt: ${event.attemptCount}"
                if( event.lastResult!=null )
                    msg += "; response: ${event.lastResult}"
                if( event.lastFailure != null )
                    msg += "; exception: [${event.lastFailure.class.name}] ${event.lastFailure.message}"
                log.debug(msg)
            }
        }
        return RetryPolicy.<T>builder()
                .handleIf(cond)
                .handleResultIf(handle)
                .withBackoff(cfg.delay.toMillis(), cfg.maxDelay.toMillis(), ChronoUnit.MILLIS)
                .withMaxAttempts(cfg.maxAttempts)
                .withJitter(cfg.jitter)
                .onRetry(listener)
                .build()
    }

    protected <T> HttpResponse<T> safeApply(CheckedSupplier action) {
        final retryOnException = (e -> e instanceof IOException) as Predicate<? extends Throwable>
        final retryOnStatusCode = ((HttpResponse<T> resp) -> resp.statusCode() in SERVER_ERRORS) as Predicate<HttpResponse<T>>
        final policy = retryPolicy(retryOnException, retryOnStatusCode)
        return Failsafe.with(policy).get(action)
    }

    static private final List<Integer> SERVER_ERRORS = [429,500,502,503,504]

    protected HttpResponse<String> httpSend(HttpRequest req)  {
        return safeApply(() -> httpClient.send(req, HttpResponse.BodyHandlers.ofString()))
    }
}
