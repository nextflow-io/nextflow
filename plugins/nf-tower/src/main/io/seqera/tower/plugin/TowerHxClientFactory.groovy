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

package io.seqera.tower.plugin

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.file.http.XAuthProvider
import nextflow.file.http.XAuthRegistry
import nextflow.trace.TraceObserverV2
import nextflow.util.Duration

import java.net.http.HttpClient

/**
 * Create and register the Tower observer instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerHxClientFactory {

    static private final String TOKEN_PREFIX = '@token:'

    @Memoized
    static HxClient httpClient(String accessToken, String refreshToken, String endpoint, TowerRetryPolicy retryPolicy) {
        assert accessToken
        final builder = HxClient.newBuilder()
        // auth settings
        setupClientAuth(builder, accessToken, refreshToken, endpoint)
        // retry settings
        builder
            .retryConfig(retryPolicy)
            .followRedirects(HttpClient.Redirect.NORMAL)
            .version(HttpClient.Version.HTTP_1_1)
            .connectTimeout(java.time.Duration.ofSeconds(60))
            .build()
    }

    private static void setupClientAuth(HxClient.Builder config, String token, String refreshToken, String endpoint) {
        // check for plain jwt token

        final refreshUrl = refreshToken ? "$endpoint/oauth/access_token" : null
        if( token.count('.')==2 ) {
            config.bearerToken(token)
            config.refreshToken(refreshToken)
            config.refreshTokenUrl(refreshUrl)
            config.refreshCookiePolicy(CookiePolicy.ACCEPT_ALL)
            return
        }

        // try checking personal access token
        try {
            final plain = new String(token.decodeBase64())
            final p = plain.indexOf('.')
            if( p!=-1 && new JsonSlurper().parseText(  plain.substring(0, p) )  ) {
                // ok this is bearer token
                config.bearerToken(token)
                // setup the refresh
                config.refreshToken(refreshToken)
                config.refreshTokenUrl(refreshUrl)
                config.refreshCookiePolicy(CookiePolicy.ACCEPT_ALL)
                return
            }
        }
        catch ( Exception e ) {
            log.trace "Enable to set bearer token ~ Reason: $e.message"
        }

        // fallback on simple token
        config.basicAuth(TOKEN_PREFIX + token)
    }
}
