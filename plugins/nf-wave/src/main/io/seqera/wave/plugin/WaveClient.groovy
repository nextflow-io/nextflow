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

package io.seqera.wave.plugin

import static io.seqera.wave.util.DockerHelper.*

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Path
import java.time.Duration
import java.time.Instant
import java.time.OffsetDateTime
import java.time.temporal.ChronoUnit
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.function.Predicate

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.util.concurrent.RateLimiter
import com.google.common.util.concurrent.UncheckedExecutionException
import com.google.gson.Gson
import com.google.gson.reflect.TypeToken
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.seqera.wave.api.BuildStatusResponse
import io.seqera.wave.api.ContainerStatus
import io.seqera.wave.api.ContainerStatusResponse
import io.seqera.wave.api.PackagesSpec
import io.seqera.wave.plugin.config.TowerConfig
import io.seqera.wave.plugin.config.WaveConfig
import io.seqera.wave.plugin.exception.BadResponseException
import io.seqera.wave.plugin.exception.UnauthorizedException
import io.seqera.wave.plugin.packer.Packer
import nextflow.Session
import nextflow.SysEnv
import nextflow.container.inspect.ContainerInspectMode
import nextflow.container.resolver.ContainerInfo
import nextflow.exception.ProcessUnrecoverableException
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

    @Canonical
    static class Handle {
        final SubmitContainerTokenResponse response
        final Instant createdAt
        int iteration
    }

    final static public String DEFAULT_S5CMD_AMD64_URL = 'https://nf-xpack.seqera.io/s5cmd/linux_amd64_2.2.2.json'
    final static public String DEFAULT_S5CMD_ARM64_URL = 'https://nf-xpack.seqera.io/s5cmd/linux_arm64_2.2.2.json'

    private static Logger log = LoggerFactory.getLogger(WaveClient)

    public static final List<String> DEFAULT_CONDA_CHANNELS = ['conda-forge','bioconda']

    private static final String DEFAULT_DOCKER_PLATFORM = 'linux/amd64'

    final private HttpClient httpClient

    final private WaveConfig config

    final private FusionConfig fusion

    final private TowerConfig tower

    final private Packer packer

    final private String endpoint

    private Cache<String, SubmitContainerTokenResponse> cache

    private Map<String,Handle> responses = new ConcurrentHashMap<>()

    private Session session

    private volatile String accessToken

    private volatile String refreshToken

    private CookieManager cookieManager

    private List<String> condaChannels

    final private String waveRegistry

    final private boolean awsFargate

    final private URL s5cmdConfigUrl

    final private RateLimiter limiter

    WaveClient(Session session) {
        this.session = session
        this.config = new WaveConfig(session.config.wave as Map ?: Collections.emptyMap(), SysEnv.get())
        this.fusion = new FusionConfig(session.config.fusion as Map ?: Collections.emptyMap(), SysEnv.get())
        this.tower = new TowerConfig(session.config.tower as Map ?: Collections.emptyMap(), SysEnv.get())
        this.awsFargate = WaveFactory.isAwsBatchFargateMode(session.config)
        this.s5cmdConfigUrl = session.config.navigate('wave.s5cmdConfigUrl') as URL
        this.endpoint = config.endpoint()
        this.condaChannels = session.getCondaConfig()?.getChannels() ?: DEFAULT_CONDA_CHANNELS
        log.debug "Wave config: $config"
        this.packer = new Packer().withPreserveTimestamp(config.preserveFileTimestamp())
        this.waveRegistry = new URI(endpoint).getAuthority()
        this.limiter = RateLimiter.create( config.httpOpts().maxRate().rate  )
        // create cache
        this.cache = CacheBuilder<String, Handle>
            .newBuilder()
            .expireAfterWrite(config.tokensCacheMaxDuration().toSeconds(), TimeUnit.SECONDS)
            .build()
        // the cookie manager
        cookieManager = new CookieManager()
        // create http client
        this.httpClient = newHttpClient()
    }

    /* only for testing */
    protected WaveClient() { }

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
        ContainerConfig containerConfig = assets.containerConfig ?: new ContainerConfig()
        // prepend the bundle layer
        if( assets.moduleResources!=null && assets.moduleResources.hasEntries() ) {
            containerConfig.prependLayer(makeLayer(assets.moduleResources))
        }
        // prepend project resources bundle
        if( assets.projectResources!=null && assets.projectResources.hasEntries() ) {
            containerConfig.prependLayer(makeLayer(assets.projectResources))
        }

        if( !assets.containerImage && !assets.containerFile && !assets.packagesSpec )
            throw new IllegalArgumentException("Wave container request requires at least a image or container file or packages spec to build")

        if( assets.containerImage && assets.containerFile )
            throw new IllegalArgumentException("Wave container image and container file cannot be specified in the same request")

        if( assets.containerImage && assets.packagesSpec )
            throw new IllegalArgumentException("Wave container image and packages spec cannot be specified in the same request")

        if( assets.containerFile && assets.packagesSpec )
            throw new IllegalArgumentException("Wave containerFile file and packages spec cannot be specified in the same request")

        if( config.mirrorMode() && config.freezeMode() )
            throw new IllegalArgumentException("Wave configuration setting 'wave.mirror' and 'wave.freeze' conflicts each other")

        if( config.mirrorMode() && !config.buildRepository() )
            throw new IllegalArgumentException("Wave configuration setting 'wave.mirror' requires the use of 'wave.build.repository' to define the target registry")

        if( config.mirrorMode() && !assets.containerImage )
            throw new IllegalArgumentException("Invalid container mirror operation - missing source container")

        if( config.mirrorMode() && containerConfig ) {
            log.warn1("Wave configuration setting 'wave.mirror' conflicts with the use of module bundles - ignoring custom config for container: $assets.containerImage")
            containerConfig = null
        }

        return new SubmitContainerTokenRequest(
                containerImage: assets.containerImage,
                containerPlatform: assets.containerPlatform,
                containerConfig: containerConfig,
                containerFile: assets.dockerFileEncoded(),
                packages: assets.packagesSpec,
                buildRepository: config().buildRepository(),
                cacheRepository: config.cacheRepository(),
                timestamp: OffsetDateTime.now().toString(),
                fingerprint: assets.fingerprint(),
                freeze: config.freezeMode(),
                format: assets.singularity ? 'sif' : null,
                dryRun: ContainerInspectMode.dryRun(),
                mirror: config.mirrorMode(),
                scanMode: config.scanMode(),
                scanLevels: config.scanAllowedLevels()
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
                dryRun: ContainerInspectMode.dryRun(),
                mirror: config.mirrorMode(),
                scanMode: config.scanMode(),
                scanLevels: config.scanAllowedLevels()
        )
        return sendRequest(request)
    }

    private void checkLimiter() {
        final ts = System.currentTimeMillis()
        try {
            limiter.acquire()
        } finally {
            final delta = System.currentTimeMillis()-ts
            if( delta>0 )
                log.debug "Request limiter blocked ${Duration.ofMillis(delta)}"
        }
    }


    SubmitContainerTokenResponse sendRequest(SubmitContainerTokenRequest request) {
        checkLimiter()
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
        final uri = URI.create("${endpoint}/v1alpha2/container")
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
                throw new BadResponseException("Wave invalid response: POST ${uri} [${resp.statusCode()}] ${resp.body()}")
        }
        catch (IOException e) {
            throw new IllegalStateException("Unable to connect Wave service: $endpoint")
        }
    }

    protected ContainerStatusResponse jsonToContainerStatusResponse(String body) {
        final obj = new JsonSlurper().parseText(body) as Map
        return new ContainerStatusResponse(
            obj.id as String,
            obj.status as ContainerStatus,
            obj.buildId as String,
            obj.mirrorId as String,
            obj.scanId as String,
            obj.vulnerabilities as Map<String,Integer>,
            obj.succeeded as Boolean,
            obj.reason as String,
            obj.detailsUri as String,
            Instant.parse(obj.creationTime as String),
            null
        )
    }

    protected BuildStatusResponse jsonToBuildStatusResponse(String body) {
        final obj = new JsonSlurper().parseText(body) as Map
        new BuildStatusResponse(
            obj.id as String,
            obj.status as BuildStatusResponse.Status,
            obj.startTime ? Instant.parse(obj.startTime as String) : null,
            obj.duration ? Duration.ofMillis(obj.duration as double * 1_000 as long) : null,
            obj.succeeded as Boolean
        )
    }

    protected SubmitContainerTokenResponse jsonToSubmitResponse(String body) {
        final type = new TypeToken<SubmitContainerTokenResponse>(){}.getType()
        return new Gson().fromJson(body, type)
    }

    protected ContainerConfig jsonToContainerConfig(String json) {
        final type = new TypeToken<ContainerConfig>(){}.getType()
        return new Gson().fromJson(json, type)
    }

    protected URL defaultFusionUrl(String platform) {
        final isArm = platform.tokenize('/')?.contains('arm64')
        return isArm
                ? new URL(FusionConfig.DEFAULT_FUSION_ARM64_URL)
                : new URL(FusionConfig.DEFAULT_FUSION_AMD64_URL)
    }

    protected URL defaultS5cmdUrl(String platform) {
        final isArm = platform.tokenize('/')?.contains('arm64')
        return isArm
            ? new URL(DEFAULT_S5CMD_ARM64_URL)
            : new URL(DEFAULT_S5CMD_AMD64_URL)
    }

    ContainerConfig resolveContainerConfig(String platform = DEFAULT_DOCKER_PLATFORM) {
        final urls = new ArrayList<URL>(config.containerConfigUrl())
        if( fusion.enabled() ) {
            final fusionUrl = fusion.containerConfigUrl() ?: defaultFusionUrl(platform)
            urls.add(fusionUrl)
        }
        if( awsFargate ) {
            final s5cmdUrl = s5cmdConfigUrl ?: defaultS5cmdUrl(platform)
            urls.add(s5cmdUrl)
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
        // get the architecture
        final arch = task.config.getArchitecture()
        final dockerArch = arch? arch.dockerArch : DEFAULT_DOCKER_PLATFORM
        // compose the request attributes
        def attrs = new HashMap<String,String>()
        attrs.container = containerImage
        attrs.conda = task.config.conda as String
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
        return resolveAssets0(attrs, bundle, singularity, dockerArch)
    }

    protected WaveAssets resolveAssets0(Map<String,String> attrs, ResourcesBundle bundle, boolean singularity, String dockerArch) {

        final scriptType = singularity ? 'singularityfile' : 'dockerfile'
        String containerScript = attrs.get(scriptType)
        final containerImage = attrs.container
        PackagesSpec packagesSpec = null

        /*
         * If 'conda' directive is specified use it to create a container file
         * to assemble the target container
         */
        if( attrs.conda ) {
            if( containerScript )
                throw new IllegalArgumentException("Unexpected conda and $scriptType conflict while resolving wave container")

            // map the recipe to a dockerfile
            if( isCondaRemoteFile(attrs.conda) ) {
                packagesSpec = new PackagesSpec()
                    .withType(PackagesSpec.Type.CONDA)
                    .withChannels(condaChannels)
                    .withCondaOpts(config.condaOpts())
                    .withEntries(List.of(attrs.conda))
            }
            else {
                if( isCondaLocalFile(attrs.conda) ) {
                    // 'conda' attribute is the path to the local conda environment
                    // note: ignore the 'channels' attribute because they are supposed to be provided by the conda file
                    final condaFile = condaFileFromPath(attrs.conda, null)
                    packagesSpec = new PackagesSpec()
                        .withType(PackagesSpec.Type.CONDA)
                        .withCondaOpts(config.condaOpts())
                        .withEnvironment(condaFile.bytes.encodeBase64().toString())
                }
                else {
                    // 'conda' attributes is resolved as the conda packages to be used
                    packagesSpec = new PackagesSpec()
                        .withType(PackagesSpec.Type.CONDA)
                        .withChannels(condaChannels)
                        .withCondaOpts(config.condaOpts())
                        .withEntries(condaPackagesToList(attrs.conda))
                }

            }
        }

        /*
         * The process should declare at least a container image name via 'container' directive
         * or a dockerfile file to build, otherwise there's no job to be done by wave
         */
        if( !containerScript && !containerImage && !packagesSpec ) {
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
                    packagesSpec,
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
            log.trace "Wave fingerprint: $key; assets: $assets"
            // get from cache or submit a new request
            final resp = cache.get(key, () -> {
                final ret = sendRequest(assets);
                responses.put(key,new Handle(ret,Instant.now()));
                return ret
            })
            return new ContainerInfo(assets.containerImage, resp.targetImage, key)
        }
        catch ( UncheckedExecutionException e ) {
            throw e.cause
        }
    }

    protected boolean checkContainerCompletion(Handle handle) {
        final long maxAwait = config.buildMaxDuration().toMillis()
        final startTime = handle.createdAt.toEpochMilli()
        final containerImage = handle.response.targetImage
        final requestId = handle.response.requestId
        final resp = containerStatus(requestId)
        if( resp.status==ContainerStatus.DONE ) {
            if( resp.succeeded )
                return true
            def msg = "Wave provisioning for container '${containerImage}' did not complete successfully"
            if( resp.reason )
                msg += "\n- Reason: ${resp.reason}"
            if( resp.detailsUri )
                msg += "\n- Find out more here: ${resp.detailsUri}"
            throw new ProcessUnrecoverableException(msg)
        }
        if( System.currentTimeMillis()-startTime > maxAwait ) {
            final msg = "Wave provisioning for container '${containerImage}' is exceeding max allowed duration (${config.buildMaxDuration()}) - check details here: ${endpoint}/view/containers/${requestId}"
            throw new ProcessUnrecoverableException(msg)
        }
        // this is expected to be invoked ~ every seconds, therefore
        // print an info message after 10 seconds or every 200 seconds
        if( ((handle.iteration++)-10) % 200 == 0 ) {
            log.info "Awaiting container provisioning: $containerImage"
        }
        return false
    }

    protected boolean checkBuildCompletion(Handle handle) {
        final long maxAwait = config.buildMaxDuration().toMillis()
        final startTime = handle.createdAt.toEpochMilli()
        final containerImage = handle.response.targetImage
        final buildId = handle.response.buildId
        final resp = buildStatus(buildId)
        if( resp.status==BuildStatusResponse.Status.COMPLETED ) {
            if( resp.succeeded )
                return true
            final msg = "Wave provisioning for container '${containerImage}' did not complete successfully - check details here: ${endpoint}/view/builds/${buildId}"
            throw new ProcessUnrecoverableException(msg)
        }
        if( System.currentTimeMillis()-startTime > maxAwait ) {
            final msg = "Wave provisioning for container '${containerImage}' is exceeding max allowed duration (${config.buildMaxDuration()}) - check details here: ${endpoint}/view/builds/${buildId}"
            throw new ProcessUnrecoverableException(msg)
        }
        // this is expected to be invoked ~ every seconds, therefore
        // print an info message after 10 seconds or every 200 seconds
        if( ((handle.iteration++)-10) % 200 == 0 ) {
            log.info "Awaiting container provisioning: $containerImage"
        }
        return false
    }

    boolean isContainerReady(String key) {
        final handle = responses.get(key)
        if( !handle )
            throw new IllegalStateException("Unable to find any container with key: $key")
        final resp = handle.response
        if( resp.requestId ) {
            return resp.succeeded
                    ? true
                    : checkContainerCompletion(handle)
        }
        if( resp.buildId && !resp.cached )
            return checkBuildCompletion(handle)
        else
            return true
    }

    protected static int randomRange(int min, int max) {
        assert min<max
        Random rand = new Random();
        return rand.nextInt((max - min) + 1) + min;
    }

    protected ContainerStatusResponse containerStatus(String requestId) {
        final String statusEndpoint = endpoint + "/v1alpha2/container/${requestId}/status";
        final HttpRequest req = HttpRequest.newBuilder()
            .uri(URI.create(statusEndpoint))
            .headers("Content-Type","application/json")
            .GET()
            .build();

        final HttpResponse<String> resp = httpSend(req);
        log.debug("Wave container status response: statusCode={}; body={}", resp.statusCode(), resp.body())
        if( resp.statusCode()==200 ) {
            return jsonToContainerStatusResponse(resp.body())
        }
        else {
            String msg = String.format("Wave invalid response: GET %s [%s] %s", statusEndpoint, resp.statusCode(), resp.body());
            throw new BadResponseException(msg)
        }
    }

    @Deprecated
    protected BuildStatusResponse buildStatus(String buildId) {
        final String statusEndpoint = endpoint + "/v1alpha1/builds/${buildId}/status";
        final HttpRequest req = HttpRequest.newBuilder()
            .uri(URI.create(statusEndpoint))
            .headers("Content-Type","application/json")
            .GET()
            .build();

        final HttpResponse<String> resp = httpSend(req);
        log.debug("Wave build status response: statusCode={}; body={}", resp.statusCode(), resp.body())
        if( resp.statusCode()==200 ) {
            return jsonToBuildStatusResponse(resp.body())
        }
        else {
            String msg = String.format("Wave invalid response: GET %s [%s] %s", statusEndpoint, resp.statusCode(), resp.body());
            throw new BadResponseException(msg)
        }
    }

    static protected boolean isCondaLocalFile(String value) {
        if( value.contains('\n') )
            return false
        if( value.startsWith('http://') || value.startsWith('https://') )
            return false
        if( value.startsWith('/') && !value.contains('\n') )
            return true
        return value.endsWith('.yaml') || value.endsWith('.yml') || value.endsWith('.txt')
    }

    static protected boolean isCondaRemoteFile(String value) {
        value.startsWith('http://') || value.startsWith('https://')
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
            void accept(ExecutionAttemptedEvent event) throws Throwable {
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
