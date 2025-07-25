package io.seqera.tower.plugin

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration
import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.concurrent.Executors
import java.util.function.Predicate

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.util.concurrent.UncheckedExecutionException
import com.google.gson.Gson
import com.google.gson.JsonSyntaxException
import dev.failsafe.Failsafe
import dev.failsafe.FailsafeException
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.exception.BadResponseException
import io.seqera.tower.plugin.exception.UnauthorizedException
import io.seqera.tower.plugin.exchange.GetLicenseTokenRequest
import io.seqera.tower.plugin.exchange.GetLicenseTokenResponse
import io.seqera.util.trace.TraceUtils
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.exception.ReportWarningException
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionToken
import nextflow.platform.PlatformHelper
import nextflow.plugin.Priority
import nextflow.serde.gson.GsonEncoder
import nextflow.util.Threads
import org.pf4j.Extension
/**
 * Environment provider for Platform-specific environment variables.
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@Slf4j
@Extension
@CompileStatic
@Priority(-10)
class TowerFusionToken implements FusionToken {

    // The path relative to the Platform endpoint where license-scoped JWT tokens are obtained
    private static final String LICENSE_TOKEN_PATH = 'license/token/'

    // Server errors that should trigger a retry
    private static final List<Integer> SERVER_ERRORS = [408, 429, 500, 502, 503, 504]

    // Default connection timeout for HTTP requests
    private static final Duration DEFAULT_CONNECTION_TIMEOUT = Duration.of(30, ChronoUnit.SECONDS)

    // Default retry policy settings for HTTP requests: delay, max delay, attempts, and jitter
    private static final Duration DEFAULT_RETRY_POLICY_DELAY = Duration.of(450, ChronoUnit.MILLIS)
    private static final Duration DEFAULT_RETRY_POLICY_MAX_DELAY = Duration.of(90, ChronoUnit.SECONDS)
    private static final int DEFAULT_RETRY_POLICY_MAX_ATTEMPTS = 10
    private static final double DEFAULT_RETRY_POLICY_JITTER = 0.5

    private CookieManager cookieManager = new CookieManager()

    // The HttpClient instance used to send requests
    private final HttpClient httpClient = newDefaultHttpClient()

    // The RetryPolicy instance used to retry requests
    private final RetryPolicy retryPolicy = newDefaultRetryPolicy(SERVER_ERRORS)

    // Time-to-live for cached tokens
    private Duration tokenTTL = Duration.of(1, ChronoUnit.HOURS)

    // Cache used for storing license tokens
    private Cache<String, GetLicenseTokenResponse> tokenCache = CacheBuilder.newBuilder()
        .expireAfterWrite(tokenTTL)
        .build()

    // Platform endpoint to use for requests
    private String endpoint

    // Platform access token to use for requests
    private volatile String accessToken

    private volatile String refreshToken

    // Platform workflowId
    private String workspaceId

    // Platform workflowId
    private String workflowId

    TowerFusionToken() {
        final config = PlatformHelper.config()
        final env = SysEnv.get()
        this.endpoint = PlatformHelper.getEndpoint(config, env)
        this.accessToken = PlatformHelper.getAccessToken(config, env)
        this.refreshToken = PlatformHelper.getRefreshToken(config, env)
        this.workflowId = env.get('TOWER_WORKFLOW_ID')
        this.workspaceId = PlatformHelper.getWorkspaceId(config, env)
    }

    protected void validateConfig() {
        if( !endpoint )
            throw new IllegalArgumentException("Missing Seqera Platform endpoint")
        if( !accessToken )
            throw new IllegalArgumentException("Missing Seqera Platform access token")
    }

    /**
     * Return any environment variables relevant to Fusion execution. This method is called
     * by {@link nextflow.fusion.FusionEnvProvider#getEnvironment} to determine which
     * environment variables are needed for the current run.
     *
     * @param scheme The scheme for which the environment variables are needed (currently unused)
     * @param config The Fusion configuration object
     * @return A map of environment variables
     */
    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        try {
            return getEnvironment0(scheme, config)
        }
        catch (Exception e) {
            final msg = "Unable to validate Fusion license - reason: ${e.message}"
            throw new ReportWarningException(msg, 'getFusionLicenseException', e)
        }
    }

    protected Map<String,String> getEnvironment0(String scheme, FusionConfig config) {
        validateConfig()
        final product = config.sku()
        final version = config.version()
        final token = getLicenseToken(product, version)
        return Map.of('FUSION_LICENSE_TOKEN', token)
    }

    /**
     * Send a request to Platform to obtain a license-scoped JWT for Fusion. The request is authenticated using the
     * Platform access token provided in the configuration of the current session.
     *
     * @throws AbortOperationException if a Platform access token cannot be found
     *
     * @return The signed JWT token
     */
    protected String getLicenseToken(String product, String version) {
        final req = new GetLicenseTokenRequest(
            product: product,
            version: version ?: 'unknown',
            workflowId: workflowId,
            workspaceId: workspaceId
        )
        final key = '${product}-${version}'
        try {
            final now = Instant.now()
            int i=0
            while( i++<2 ) {
                final resp = tokenCache.get(key, () -> sendRequest(req))
                if( resp.error )
                    throw resp.error
                // Check if the cached response has expired
                // It's needed because the JWT token TTL in the cache (1 hour) and its expiration date (e.g. 1 day?) are not sync'ed,
                // so it could happen that we get a token from the cache which was valid at the time of insertion but is now expired.
                if( resp.expiresAt.isBefore(now) ) {
                    log.debug "Cached token already expired; refreshing"
                    tokenCache.invalidate(key)
                }
                else
                    return resp.signedToken
            }
        }
        catch (UncheckedExecutionException e) {
            // most likely the exception is thrown for the lack of license
            // to avoid to keep requesting it, and error response is added to the cache
            tokenCache.put(key, new GetLicenseTokenResponse(error: e.cause))
            throw e.cause
        }
    }

    /**************************************************************************
     * Helper methods
     *************************************************************************/

    /**
     * Create a new HttpClient instance with default settings
     * @return The new HttpClient instance
     */
    private HttpClient newDefaultHttpClient() {
        final builder = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_1_1)
            .followRedirects(HttpClient.Redirect.NEVER)
            .cookieHandler(cookieManager)
            .connectTimeout(DEFAULT_CONNECTION_TIMEOUT)
        // use virtual threads executor if enabled
        if ( Threads.useVirtual() ) {
            builder.executor(Executors.newVirtualThreadPerTaskExecutor())
        }
        // build and return the new client
        return builder.build()
    }

    /**
     * Create a new RetryPolicy instance with default settings and the given list of retryable errors. With this policy,
     * a request is retried on IOExceptions and any server errors defined in errorsToRetry. The number of retries, delay,
     * max delay, and jitter are controlled by the corresponding values defined at class level.
     *
     * @return The new RetryPolicy instance
     */
    private static <T> RetryPolicy<HttpResponse<T>> newDefaultRetryPolicy(List<Integer> errorsToRetry) {

        final retryOnException = (e -> e instanceof IOException) as Predicate<? extends Throwable>
        final retryOnStatusCode = ((HttpResponse<T> resp) -> resp.statusCode() in errorsToRetry) as Predicate<HttpResponse<T>>

        final listener = new EventListener<ExecutionAttemptedEvent<HttpResponse<T>>>() {
            @Override
            void accept(ExecutionAttemptedEvent event) throws Throwable {
                def msg = "connection failure - attempt: ${event.attemptCount}"
                if (event.lastResult != null)
                    msg += "; response: ${event.lastResult}"
                if (event.lastFailure != null)
                    msg += "; exception: [${event.lastFailure.class.name}] ${event.lastFailure.message}"
                log.debug(msg)
            }
        }
        return RetryPolicy.<HttpResponse<T>> builder()
            .handleIf(retryOnException)
            .handleResultIf(retryOnStatusCode)
            .withBackoff(DEFAULT_RETRY_POLICY_DELAY.toMillis(), DEFAULT_RETRY_POLICY_MAX_DELAY.toMillis(), ChronoUnit.MILLIS)
            .withMaxAttempts(DEFAULT_RETRY_POLICY_MAX_ATTEMPTS)
            .withJitter(DEFAULT_RETRY_POLICY_JITTER)
            .onRetry(listener)
            .build()
    }

    /**
     * Send an HTTP request and return the response. This method automatically retries the request according to the
     * given RetryPolicy.
     *
     * @param req The HttpRequest to send
     * @return The HttpResponse received
     */
    private <T> HttpResponse<String> safeHttpSend(HttpRequest req) {
        try {
            safeApply(req)
        }
        catch (FailsafeException e) {
            throw e.cause
        }
    }

    private <T> HttpResponse<String> safeApply(HttpRequest req) {
        return Failsafe.with(retryPolicy).get(
            () -> {
                log.debug "Http request: method=${req.method()}; uri=${req.uri()}; request=${req}"
                final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
                log.debug "Http response: statusCode=${resp.statusCode()}; body=${resp.body()}"
                return resp
            } as CheckedSupplier
        ) as HttpResponse<String>
    }

    /**
     * Create a {@link HttpRequest} representing a {@link GetLicenseTokenRequest} object
     *
     * @param req The LicenseTokenRequest object
     * @return The resulting HttpRequest object
     */
    private HttpRequest makeHttpRequest(GetLicenseTokenRequest req) {
        final gson = new GsonEncoder<GetLicenseTokenRequest>() {}
        final body = HttpRequest.BodyPublishers.ofString( gson.encode(req) )
        return HttpRequest.newBuilder()
            .uri(URI.create("${endpoint}/${LICENSE_TOKEN_PATH}").normalize())
            .header('Content-Type', 'application/json')
            .header('Traceparent', TraceUtils.rndTrace())
            .header('Authorization', "Bearer ${accessToken}")
            .POST(body)
            .build()
    }

    /**
     * Serialize a {@link GetLicenseTokenRequest} object into a JSON string
     *
     * @param req The LicenseTokenRequest object
     * @return The resulting JSON string
     */
    private static String serializeToJson(GetLicenseTokenRequest req) {
        return new Gson().toJson(req)
    }

    /**
     * Parse a JSON string into a {@link GetLicenseTokenResponse} object
     *
     * @param json The String containing the JSON representation of the LicenseTokenResponse object
     * @return The resulting LicenseTokenResponse object
     *
     * @throws JsonSyntaxException if the JSON string is not well-formed
     */
    protected static GetLicenseTokenResponse parseLicenseTokenResponse(String json) throws JsonSyntaxException {
        final gson = new GsonEncoder<GetLicenseTokenResponse>() {}
        return gson.decode(json)
    }

    /**
     * Request a license token from Platform.
     *
     * @param request The LicenseTokenRequest object
     * @return The LicenseTokenResponse object
     */
    private GetLicenseTokenResponse sendRequest(GetLicenseTokenRequest request) {
        return sendRequest0(request, 1)
    }

    private GetLicenseTokenResponse sendRequest0(GetLicenseTokenRequest request, int attempt) {

        final httpReq = makeHttpRequest(request)

        try {
            final resp = safeHttpSend(httpReq)

            if( resp.statusCode() == 200 ) {
                final ret = parseLicenseTokenResponse(resp.body())
                return ret
            }

            if( resp.statusCode() == 401 ) {
                final shouldRetry = accessToken
                    && refreshToken
                    && attempt==1
                    && refreshJwtToken0(refreshToken)
                if( shouldRetry ) {
                    return sendRequest0(request, attempt+1)
                }
                else
                    throw new UnauthorizedException("Unauthorized [401] - Verify you have provided a Seqera Platform valid access token")
            }

            throw new BadResponseException("Invalid response: ${httpReq.method()} ${httpReq.uri()} [${resp.statusCode()}] ${resp.body()}")
        }
        catch (IOException e) {
            throw new IllegalStateException("Unable to send request to '${httpReq.uri()}' : ${e.message}")
        }
    }

    protected boolean refreshJwtToken0(String refresh) {
        log.debug "Token refresh request >> $refresh"

        final req = HttpRequest.newBuilder()
            .uri(new URI("${endpoint}/oauth/access_token"))
            .headers('Content-Type',"application/x-www-form-urlencoded")
            .POST(HttpRequest.BodyPublishers.ofString("grant_type=refresh_token&refresh_token=${URLEncoder.encode(refresh, 'UTF-8')}"))
            .build()

        final resp = safeHttpSend(req)
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
}
