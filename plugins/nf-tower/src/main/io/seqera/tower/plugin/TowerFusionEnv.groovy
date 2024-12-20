package io.seqera.tower.plugin

import com.google.gson.Gson
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.event.EventListener
import dev.failsafe.event.ExecutionAttemptedEvent
import dev.failsafe.function.CheckedSupplier
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.exception.BadResponseException
import io.seqera.tower.plugin.exception.UnauthorizedException
import io.seqera.tower.plugin.exchange.LicenseTokenRequest
import io.seqera.tower.plugin.exchange.LicenseTokenResponse
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionEnv
import nextflow.platform.PlatformHelper
import nextflow.util.Threads
import org.pf4j.Extension

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration
import java.time.temporal.ChronoUnit
import java.util.concurrent.Executors
import java.util.function.Predicate

/**
 * Environment provider for Platform-specific environment variables.
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@Slf4j
@Extension
@CompileStatic
class TowerFusionEnv implements FusionEnv {

    // The path relative to the Platform endpoint where license-scoped JWT tokens are obtained
    private static final String LICENSE_TOKEN_PATH = 'license/token/'

    // Server errors that should trigger a retry
    private static final List<Integer> SERVER_ERRORS = [429, 500, 502, 503, 504]

    // Default connection timeout for HTTP requests
    private static final Duration DEFAULT_CONNECTION_TIMEOUT = Duration.of(30, ChronoUnit.SECONDS)

    // Default retry policy settings for HTTP requests: delay, max delay, attempts, and jitter
    private static final Duration DEFAULT_RETRY_POLICY_DELAY = Duration.of(450, ChronoUnit.MILLIS)
    private static final Duration DEFAULT_RETRY_POLICY_MAX_DELAY = Duration.of(90, ChronoUnit.SECONDS)
    private static final int DEFAULT_RETRY_POLICY_MAX_ATTEMPTS = 10
    private static final double DEFAULT_RETRY_POLICY_JITTER = 0.5

    // The HttpClient instance used to send requests
    private final HttpClient httpClient = newDefaultHttpClient()

    // The RetryPolicy instance used to retry requests
    private final RetryPolicy retryPolicy = newDefaultRetryPolicy(SERVER_ERRORS)

    // Nextflow session
    private final Session session

    // Platform endpoint to use for requests
    private final String endpoint

    // Platform access token to use for requests
    private final String accessToken

    /**
     * Constructor for the class. It initializes the session, endpoint, and access token.
     */
    TowerFusionEnv() {
        this.session = Global.session as Session
        final towerConfig = session.config.navigate('tower') as Map ?: [:]
        final env = SysEnv.get()
        this.endpoint = PlatformHelper.getEndpoint(towerConfig, env)
        this.accessToken = PlatformHelper.getAccessToken(towerConfig, env)
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

        final product = config.sku()
        final version = config.version()

        try {
            final token = getLicenseToken(product, version)
            return [
                FUSION_LICENSE_TOKEN: token,
            ]
        } catch (Exception e) {
            log.warn("Error retrieving Fusion license information: ${e.message}")
            return Map.of()
        }
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
        // FIXME(amiranda): Find out how to obtain the product and version
        // Candidate: FusionConfig?

        if (accessToken == null) {
            throw new AbortOperationException("Missing personal access token -- Make sure there's a variable TOWER_ACCESS_TOKEN in your environment")
        }

        final req = HttpRequest.newBuilder()
            .uri(URI.create("${endpoint}/${LICENSE_TOKEN_PATH}").normalize())
            .header('Content-Type', 'application/json')
            .header('Authorization', "Bearer ${accessToken}")
            .POST(
                HttpRequest.BodyPublishers.ofString(
                    makeLicenseTokenRequest(product, version)
                )
            )
            .build()

        try {
            final resp = safeHttpSend(req, retryPolicy)

            if (resp.statusCode() == 200) {
                final ret = parseLicenseTokenResponse(resp)
                return ret.signedToken
            }

            if (resp.statusCode() == 401) {
                throw new UnauthorizedException("Unauthorized [401] - Verify you have provided a valid access token")
            }

            throw new BadResponseException("Invalid response: ${req.method()} ${req.uri()} [${resp.statusCode()}] ${resp.body()}")

        } catch (IOException e) {
            throw new IllegalStateException("Unable to send request to '${req.uri()}' : ${e.message}")
        }
    }

    /**************************************************************************
     * Helper methods
     *************************************************************************/

    /**
     * Create a new HttpClient instance with default settings
     * @return The new HttpClient instance
     */
    private static HttpClient newDefaultHttpClient() {
        final builder = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_1_1)
            .followRedirects(HttpClient.Redirect.NEVER)
            .cookieHandler(new CookieManager())
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
    private <T> HttpResponse<String> safeHttpSend(HttpRequest req, RetryPolicy<T> policy) {
        return Failsafe.with(policy).get(
            () -> {
                log.debug "Request: method:=${req.method()}; uri:=${req.uri()}; request:=${req}"
                final resp = httpClient.send(req, HttpResponse.BodyHandlers.ofString())
                log.debug "Response: statusCode:=${resp.statusCode()}; body:=${resp.body()}"
                return resp
            } as CheckedSupplier
        ) as HttpResponse<String>
    }

    /**
     * Create a JSON string representing a {@link LicenseTokenRequest} object
     *
     * @param product The product SKU
     * @param version The version
     * @return The resulting JSON string
     */
    private static String makeLicenseTokenRequest(String product, String version) {
        return new Gson().toJson(
            new LicenseTokenRequest(
                product: product,
                version: version
            ),
            LicenseTokenRequest.class
        )
    }

    /**
     * Parse a JSON string into a {@link LicenseTokenResponse} object
     *
     * @param stringHttpResponse The HttpResponse containing the JSON string
     * @return The resulting LicenseTokenResponse object
     */
    private static LicenseTokenResponse parseLicenseTokenResponse(HttpResponse<String> resp) {
        return new Gson().fromJson(resp.body(), LicenseTokenResponse.class)
    }
}
