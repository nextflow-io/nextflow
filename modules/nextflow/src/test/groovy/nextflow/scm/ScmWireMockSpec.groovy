/*
 * Copyright 2013-2026, Seqera Labs
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
 */

package nextflow.scm

import static com.github.tomakehurst.wiremock.client.WireMock.*
import static com.github.tomakehurst.wiremock.core.WireMockConfiguration.wireMockConfig

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.MappingBuilder
import com.github.tomakehurst.wiremock.client.ResponseDefinitionBuilder
import nextflow.SysEnv
import spock.lang.Specification

/**
 * Base for the offline SCM provider specs that run against a WireMock fake of the provider API.
 *
 * Stubs match the auth header returned by {@link #expectedAuth()} by default, so a request
 * that drops the header gets no stub match and fails the test. Pass {@link #NO_AUTH} to
 * require that no auth header is sent.
 */
abstract class ScmWireMockSpec extends Specification {

    static final String NO_AUTH = null

    static final byte[] BINARY = [0, 1, 2, -1, -2, 127] as byte[]

    WireMockServer wireMock
    String url

    def setup() {
        // isolate from SCM tokens in the host environment
        SysEnv.push([:])
        wireMock = new WireMockServer(wireMockConfig().dynamicPort())
        wireMock.start()
        url = "http://localhost:${wireMock.port()}"
    }

    def cleanup() {
        wireMock?.stop()
        SysEnv.pop()
    }

    /**
     * The auth header value sent with the default credentials of {@code provider()}
     */
    protected abstract String expectedAuth()

    /**
     * The name of the header that carries the credentials
     */
    protected String authHeader() { 'Authorization' }

    protected MappingBuilder withAuth(MappingBuilder builder, String auth=expectedAuth()) {
        return auth != null
            ? builder.withHeader(authHeader(), equalTo(auth))
            : builder.withHeader(authHeader(), absent())
    }

    protected void stub(String path, ResponseDefinitionBuilder response, String auth=expectedAuth()) {
        wireMock.stubFor(withAuth(get(urlEqualTo(path)), auth).willReturn(response))
    }

    protected void stubJson(String path, String body, String auth=expectedAuth()) {
        stub(path, okJson(body), auth)
    }

    protected void stubStatus(String path, int status, String body, String auth=expectedAuth()) {
        stub(path, aResponse().withStatus(status).withBody(body), auth)
    }
}
