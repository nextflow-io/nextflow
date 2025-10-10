package nextflow.platform

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv

/**
 * Helper methods for Platform-related operations
 *
 * @author Alberto Miranda <alberto.miranda@seqera.io>
 */
@CompileStatic
class PlatformHelper {

    /**
     * Get the configured Platform API endpoint: if the endpoint is not provided in the configuration, we fallback to the
     * environment variable `TOWER_API_ENDPOINT`. If neither is provided, we fallback to the default endpoint.
     *
     * @param opts the configuration options for Platform (e.g. `session.config.navigate('tower')`)
     * @param env the applicable environment variables
     * @return the Platform API endpoint
     */
    static String getEndpoint(Map opts, Map<String,String> env) {
        def result = opts.endpoint as String
        if( !result || result=='-' )
            result = env.get('TOWER_API_ENDPOINT') ?: 'https://api.cloud.seqera.io'
        return result.stripEnd('/')
    }

    /**
     * Get the Auth0 domain for a given Platform API endpoint
     *
     * @param endpoint the Platform API endpoint
     * @return the Auth0 domain, or null if not a cloud endpoint
     */
    static String getAuthDomain(String endpoint) {
        switch(endpoint) {
            case SysEnv.get('TOWER_AUTH_DOMAIN'):
                return SysEnv.get('TOWER_AUTH_DOMAIN')
            case 'https://api.cloud.dev-seqera.io':
                return 'seqera-development.eu.auth0.com'
            case 'https://api.cloud.stage-seqera.io':
                return 'seqera-stage.eu.auth0.com'
            case 'https://api.cloud.seqera.io':
                return 'seqera.eu.auth0.com'
            default:
                return null
        }
    }

    /**
     * Get the Auth0 client ID for a given Platform API endpoint
     *
     * @param endpoint the Platform API endpoint
     * @return the Auth0 client ID, or null if not a cloud endpoint
     */
    static String getAuthClientId(String endpoint) {
        switch(endpoint) {
            case SysEnv.get('TOWER_AUTH_ID'):
                return SysEnv.get('TOWER_AUTH_ID')
            case 'https://api.cloud.dev-seqera.io':
                return 'Ep2LhYiYmuV9hhz0dH6dbXVq0S7s7SWZ'
            case 'https://api.cloud.stage-seqera.io':
                return '60cPDjI6YhoTPjyMTIBjGtxatSUwWswB'
            case 'https://api.cloud.seqera.io':
                return 'FxCM8EJ76nNeHUDidSHkZfT8VtsrhHeL'
            default:
                return null
        }
    }

    /**
     * Return the configured Platform access token: if `TOWER_WORKFLOW_ID` is provided in the environment, it means
     * we are running in a Platform-made run and we should ONLY retrieve the token from the environment. Otherwise,
     * check the configuration or fallback to the environment. If no token is found, return null.
     *
     * @param opts the configuration options for Platform (e.g. `session.config.navigate('tower')`)
     * @param env the applicable environment variables
     * @return the Platform access token
     */
    static String getAccessToken(Map opts, Map<String,String> env) {
        final token = env.get('TOWER_WORKFLOW_ID')
            ? env.get('TOWER_ACCESS_TOKEN')
            : opts.containsKey('accessToken') ? opts.accessToken as String : env.get('TOWER_ACCESS_TOKEN')
        return token
    }

    /**
     * Return the configured Platform refresh token: if `TOWER_WORKFLOW_ID` is provided in the environment, it means
     * we are running in a Platform-made run and we should ONLY retrieve the token from the environment. Otherwise,
     * check the configuration or fallback to the environment. If no token is found, return null.
     *
     * @param opts the configuration options for Platform (e.g. `session.config.navigate('tower')`)
     * @param env the applicable environment variables
     * @return the Platform refresh token
     */
    static String getRefreshToken(Map opts, Map<String,String> env) {
        final token = env.get('TOWER_WORKFLOW_ID')
            ? env.get('TOWER_REFRESH_TOKEN')
            : opts.containsKey('refreshToken') ? opts.refreshToken as String : env.get('TOWER_REFRESH_TOKEN')
        return token
    }

    /**
     * Return the Platform Workspace ID: if `TOWER_WORKFLOW_ID` is provided in the environment, it means we are running
     * in a Platform-made run and we should ONLY retrieve the workspace ID from the environment. Otherwise, check the
     * configuration or fallback to the environment. If no workspace ID is found, return null.
     * @param opts
     * @param env
     * @return
     */
    static String getWorkspaceId(Map opts, Map<String,String> env) {
        final workspaceId = env.get('TOWER_WORKFLOW_ID')
            ? env.get('TOWER_WORKSPACE_ID')
            : opts.workspaceId as Long ?: env.get('TOWER_WORKSPACE_ID') as Long
        return workspaceId
    }

    static Map<String,Object> config() {
        session().config.navigate('tower') as Map<String,Object> ?: Map.<String,Object>of()
    }

    static private Session session() {
        Global.session as Session
    }
}
