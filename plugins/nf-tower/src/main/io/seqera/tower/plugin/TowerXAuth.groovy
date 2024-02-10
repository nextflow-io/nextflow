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

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.http.XAuthProvider

/**
 * Implements Tower authentication strategy for resources accessed
 * via {@link nextflow.file.http.XFileSystemProvider}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class TowerXAuth implements XAuthProvider {

    final private Pattern pattern
    final private String endpoint
    private String accessToken
    private String refreshToken
    private CookieManager cookieManager
    final private HttpClient httpClient

    TowerXAuth(String endpoint, String accessToken, String refreshToken) {
        this.endpoint = endpoint
        this.pattern = ~/(?i)^$endpoint\/.*$/
        this.accessToken = accessToken
        this.refreshToken = refreshToken
        //
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

    @Override
    boolean authorize(URLConnection conn) {
        final req = conn.getURL().toString()
        if( pattern.matcher(req).matches() && !conn.getRequestProperty('Authorization') ) {
            log.trace "Authorizing request connection to: $req"
            conn.setRequestProperty('Authorization', "Bearer $accessToken")
            return true
        }
        return false
    }

    boolean refreshToken(URLConnection conn) {
        if( !refreshToken || !pattern.matcher(conn.getURL().toString()).matches() ) {
            return false
        }

        final req = HttpRequest.newBuilder()
                .uri(new URI("${endpoint}/oauth/access_token"))
                .headers('Content-Type',"application/x-www-form-urlencoded")
                .POST(HttpRequest.BodyPublishers.ofString("grant_type=refresh_token&refresh_token=${URLEncoder.encode(refreshToken, 'UTF-8')}"))
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
