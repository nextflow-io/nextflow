package io.seqera.tower.plugin

import java.time.Duration
import java.time.Instant
import java.time.temporal.ChronoUnit

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.util.concurrent.UncheckedExecutionException
import com.google.gson.JsonSyntaxException
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.exception.BadResponseException
import io.seqera.tower.plugin.exception.UnauthorizedException
import io.seqera.tower.plugin.exchange.GetLicenseTokenRequest
import io.seqera.tower.plugin.exchange.GetLicenseTokenResponse
import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.exception.ReportWarningException
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionToken
import nextflow.platform.PlatformHelper
import nextflow.plugin.Priority
import nextflow.serde.gson.GsonEncoder
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

    // The TowerClient instance used to send requests
    private TowerClient client

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
        this.client = TowerFactory.client()
    }

    protected void validateConfig() {
        if( !endpoint )
            throw new IllegalArgumentException("Missing Seqera Platform endpoint")
        if( !accessToken )
            throw new IllegalArgumentException("Seqera Platform access token is required to use Fusion -- see https://docs.seqera.io/fusion/licensing for more information")
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
        final url = "${client.getEndpoint()}/${LICENSE_TOKEN_PATH}"
        final resp = client.sendHttpMessage(url, request.toMap())

        if( resp.code == 200 ) {
            final ret = parseLicenseTokenResponse(resp.message)
            return ret
        }

        if( resp.code == 401 ) {
            throw new UnauthorizedException("Unauthorized [401] - Verify you have provided a Seqera Platform valid access token")
        }

        throw new BadResponseException("Invalid response: ${url} [${resp.code}] ${resp.message} -- ${resp.cause}")

    }

}
