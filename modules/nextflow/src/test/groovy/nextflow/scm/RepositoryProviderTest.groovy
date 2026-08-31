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

import java.net.http.HttpClient
import java.net.http.HttpConnectTimeoutException
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.net.http.HttpTimeoutException
import java.nio.channels.UnresolvedAddressException
import java.time.Duration
import javax.net.ssl.SSLSession

import com.sun.net.httpserver.HttpHandler
import com.sun.net.httpserver.HttpServer
import nextflow.SysEnv
import nextflow.exception.HttpResponseLengthExceedException
import nextflow.util.HttpClientOpts
import nextflow.util.RetryConfig
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RepositoryProviderTest extends Specification {

    def 'should create repository provider object' () {

        def provider

        when:
        provider = RepositoryFactory.newRepositoryProvider(new ProviderConfig('github'),'project/x')
        then:
        provider instanceof GithubRepositoryProvider
        provider.endpointUrl == 'https://api.github.com/repos/project/x'

        when:
        provider = RepositoryFactory.newRepositoryProvider(new ProviderConfig('gitlab'),'project/y')
        then:
        provider instanceof GitlabRepositoryProvider
        provider.endpointUrl == 'https://gitlab.com/api/v4/projects/project%2Fy'

        when:
        provider = RepositoryFactory.newRepositoryProvider(new ProviderConfig('bitbucket'),'project/z')
        then:
        provider instanceof BitbucketRepositoryProvider
        provider.endpointUrl == 'https://api.bitbucket.org/2.0/repositories/project/z'

        when:
        provider = RepositoryFactory.newRepositoryProvider(new ProviderConfig('local', [path:'/user/data']),'local/w')
        then:
        provider.endpointUrl == 'file:/user/data/w'
    }

    def 'should set credentials' () {

        given:
        def config = Mock(ProviderConfig)
        def provider = Spy(RepositoryProvider)
        provider.@config = config

        when:
        provider.setCredentials('pditommaso', 'secret1')
        then:
        1 * config.setUser('pditommaso')
        1 * config.setPassword('secret1')

    }

    def 'should hide creds' () {
        given:
        def provider = Spy(RepositoryProvider)

        when:
        def result = provider.getAuthObfuscated()
        then:
        result == '-:-'

        when:
        result = provider.getAuthObfuscated()
        then:
        provider.getUser() >> 'foo1234567890'
        provider.getPassword() >> 'bar4567890'
        and:
        result == 'foo****:bar****'

    }

    def 'should auth using credentials' () {
        given:
        def provider = Spy(RepositoryProvider)
        and:
        def conn = Mock(HttpURLConnection)

        when:
        def headers = provider.getAuth()
        then:
        1 * provider.getUser() >> null
        1 * provider.hasCredentials()
        0 * conn.setRequestProperty('Authorization', _)
        and:
        !headers

        when:
        headers = provider.getAuth()
        then:
        _ * provider.getUser() >> 'foo'
        _ * provider.getPassword() >> 'bar'
        1 * provider.hasCredentials()
        and:
        headers == new String[] { 'Authorization', "Basic ${'foo:bar'.bytes.encodeBase64()}" }
    }


    def 'should validate is retryable' () {
        given:
        def provider = Spy(RepositoryProvider)
        expect:
        provider.isRetryable(ERR) == EXPECTED

        where:
        EXPECTED    | ERR
        false       | new RuntimeException()
        false       | new IOException()
        true        | new SocketException()
        false       | new SocketException(new UnresolvedAddressException())
        false       | new SocketTimeoutException()
    }

    def 'should retry connect timeouts via the http client retry condition' () {
        given:
        def provider = Spy(RepositoryProvider)

        when:
        def condition = provider.retryConfig0().getRetryCondition()
        then:
        // a connect timeout is retried, either direct or as the cause of another error
        condition.test(new HttpConnectTimeoutException('connect timed out'))
        condition.test(new IOException(new HttpConnectTimeoutException('connect timed out')))
        // a read timeout is not retried
        !condition.test(new HttpTimeoutException('request timed out'))
        // an exception without a cause does not throw a NPE
        !condition.test(new IOException('generic failure'))
    }

    def 'should not validate when max length is not configured via SysEnv' () {
        given:
        def provider = Spy(RepositoryProvider)
        def response = createMockResponseWithContentLength(1500)

        when:
        SysEnv.push([:])
        provider.checkMaxLength(response)
        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
    }

    def 'should validate and pass when content is within limit via SysEnv' () {
        given:
        def provider = Spy(RepositoryProvider)
        def response = createMockResponseWithContentLength(1500)

        when:
        SysEnv.push(['NXF_GIT_RESPONSE_MAX_LENGTH': '2000'])
        provider.checkMaxLength(response)
        then:
        noExceptionThrown()

        cleanup:
        SysEnv.pop()
    }

    def 'should validate and fail when content exceeds limit via SysEnv' () {
        given:
        def provider = Spy(RepositoryProvider)
        def response = createMockResponseWithContentLength(1500)

        when:
        SysEnv.push(['NXF_GIT_RESPONSE_MAX_LENGTH': '1000'])
        provider.checkMaxLength(response)
        then:
        thrown(HttpResponseLengthExceedException)

        cleanup:
        SysEnv.pop()
    }

    def 'should fail when the timeout value is not a duration' () {
        given:
        SysEnv.push([NXF_SCM_HTTPCLIENT_CONNECT_TIMEOUT: 'foo'])
        def provider = Spy(RepositoryProvider)

        when:
        provider.getHttpClientOpts()
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('foo')

        cleanup:
        SysEnv.pop()
    }

    def 'should prefer the http client opts set explicitly over the env variables' () {
        given:
        SysEnv.push([NXF_SCM_HTTPCLIENT_CONNECT_TIMEOUT: '5s', NXF_SCM_HTTPCLIENT_REQUEST_TIMEOUT: '250ms'])
        def provider = Spy(RepositoryProvider)
        and:
        provider.setHttpClientOpts(new HttpClientOpts(
            nextflow.util.Duration.of('30s'),
            nextflow.util.Duration.of('90s')))

        expect:
        provider.httpClientOpts.connectTimeout() == Duration.ofSeconds(30)
        provider.httpClientOpts.requestTimeout() == Duration.ofSeconds(90)

        cleanup:
        SysEnv.pop()
    }

    def 'should apply the request timeout to the http request as connect plus request' () {
        given:
        SysEnv.push(env)
        def provider = Spy(RepositoryProvider)

        when:
        def request = provider.newRequest('http://foo.com/bar')
        then:
        // HttpRequest.timeout() spans the connect phase too, so it carries the sum;
        // an unbounded request timeout applies none
        request.timeout().orElse(null) == expected

        cleanup:
        SysEnv.pop()

        where:
        env                                                                                     | expected
        [NXF_SCM_HTTPCLIENT_CONNECT_TIMEOUT: '5s', NXF_SCM_HTTPCLIENT_REQUEST_TIMEOUT: '250ms']  | Duration.ofMillis(5_250)
        [NXF_SCM_HTTPCLIENT_REQUEST_TIMEOUT: '0s']                                               | null
    }

    def 'should fail with a timeout error when the server does not reply' () {
        given:
        def server = HttpServer.create(new InetSocketAddress(0), 0)
        server.createContext('/stall', { exchange ->
            sleep 5_000
            exchange.sendResponseHeaders(200, 0)
            exchange.close()
        } as HttpHandler)
        server.start()
        and:
        SysEnv.push([NXF_SCM_HTTPCLIENT_CONNECT_TIMEOUT: '1s', NXF_SCM_HTTPCLIENT_REQUEST_TIMEOUT: '500ms'])
        def provider = Spy(RepositoryProvider)

        when:
        provider.invokeBytes("http://localhost:${server.address.port}/stall")
        then:
        def e = thrown(IOException)
        e.message.startsWith("No response received from 'http://localhost:${server.address.port}/stall'")
        e.message.contains('NXF_SCM_HTTPCLIENT_REQUEST_TIMEOUT')
        e.cause instanceof HttpTimeoutException

        cleanup:
        SysEnv.pop()
        server.stop(0)
    }

    def 'should build the http client opts from the scm config map' () {
        given:
        SysEnv.push([NXF_SCM_HTTPCLIENT_CONNECT_TIMEOUT: '5s'])

        when: 'the scm file carries an httpClient block'
        def opts = RepositoryProvider.httpClientOpts([httpClient: [connectTimeout: '30s', requestTimeout: '90s']])
        then: 'it wins over the env variable'
        opts.connectTimeout() == Duration.ofSeconds(30)
        opts.requestTimeout() == Duration.ofSeconds(90)

        when: 'the scm file has no httpClient block'
        opts = RepositoryProvider.httpClientOpts([providers: [github: [user: 'foo']]])
        then: 'the env variable applies, then the default'
        opts.connectTimeout() == Duration.ofSeconds(5)
        opts.requestTimeout() == Duration.ofSeconds(60)

        when: 'there is no scm config at all'
        opts = RepositoryProvider.httpClientOpts(null)
        then:
        opts.connectTimeout() == Duration.ofSeconds(5)

        cleanup:
        SysEnv.pop()
    }

    private createMockResponseWithContentLength(long contentLength) {
        return new HttpResponse<byte[]>() {
            @Override
            int statusCode() { return 200 }

            @Override
            HttpRequest request() { return null }

            @Override
            Optional<HttpResponse<byte[]>> previousResponse() { return Optional.empty() }

            @Override
            java.net.http.HttpHeaders headers() {
                return java.net.http.HttpHeaders.of(
                    ['Content-Length': [contentLength.toString()]],
                    (a, b) -> true
                )
            }

            @Override
            byte[] body() { return new byte[0] }

            @Override
            Optional<SSLSession> sslSession() { return Optional.empty() }

            @Override
            URI uri() { return new URI('https://api.github.com/repos/test/repo') }

            @Override
            HttpClient.Version version() { return HttpClient.Version.HTTP_1_1 }
        }
    }

    // ====== Path normalization helper method tests ======

    @Unroll
    def 'normalizePath should handle #description'() {
        expect:
        RepositoryProvider.normalizePath(INPUT) == EXPECTED

        where:
        INPUT           | EXPECTED      | description
        null            | ""            | "null input"
        ""              | ""            | "empty string"
        "/"             | ""            | "root directory slash"
        "/docs"         | "docs"        | "absolute path"
        "docs"          | "docs"        | "relative path"
        "/docs/guide"   | "docs/guide"  | "nested absolute path"
        "docs/guide"    | "docs/guide"  | "nested relative path"
        "//"            | "/"           | "double slash"
        "///docs"       | "//docs"      | "multiple leading slashes"
    }

    @Unroll
    def 'ensureAbsolutePath should handle #description'() {
        expect:
        RepositoryProvider.ensureAbsolutePath(INPUT) == EXPECTED

        where:
        INPUT           | EXPECTED        | description
        null            | "/"             | "null input"
        ""              | "/"             | "empty string"
        "/"             | "/"             | "root directory"
        "/docs"         | "/docs"         | "already absolute path"
        "docs"          | "/docs"         | "relative path"
        "/docs/guide"   | "/docs/guide"   | "nested absolute path"
        "docs/guide"    | "/docs/guide"   | "nested relative path"
        "main.nf"       | "/main.nf"      | "simple filename"
    }

    @Unroll
    def 'shouldIncludeAtDepth should handle depth=#depth basePath=#basePath entryPath=#entryPath'() {
        expect:
        RepositoryProvider.shouldIncludeAtDepth(entryPath, basePath, depth) == expected

        where:
        entryPath           | basePath    | depth | expected | description
        // Root directory tests (basePath = null, "", or "/")
        "main.nf"           | null        | 0     | true     | "immediate child in root with depth 0"
        "docs/guide.md"     | null        | 0     | false    | "nested file in root with depth 0"
        "docs/guide.md"     | null        | 1     | true     | "nested file in root with depth 1"
        "docs/sub/file.md"  | null        | 1     | false    | "deeply nested file with depth 1"
        "docs/sub/file.md"  | null        | 2     | true     | "deeply nested file with depth 2"
        "main.nf"           | ""          | 0     | true     | "immediate child with empty basePath"
        "main.nf"           | "/"         | 0     | true     | "immediate child with root basePath"

        // Subdirectory tests
        "docs/guide.md"     | "docs"      | 0     | true     | "immediate child in subdirectory"
        "docs/sub/file.md"  | "docs"      | 0     | false    | "nested file in subdirectory with depth 0"
        "docs/sub/file.md"  | "docs"      | 1     | true     | "nested file in subdirectory with depth 1"
        "docs/guide.md"     | "/docs"     | 0     | true     | "immediate child with absolute basePath"

        // Edge cases
        "docs"              | "docs"      | 0     | false    | "base directory itself should be excluded"
        "other/file.md"     | "docs"      | 0     | false    | "file outside basePath should be excluded"
        "main.nf"           | null        | -1    | true     | "unlimited depth should include everything"
        "docs/sub/deep.md"  | null        | -1    | true     | "unlimited depth with nested file"
        ""                  | null        | 0     | false    | "empty entryPath should be excluded"

        // Complex path tests
        "docs/api/index.md" | "docs"      | 1     | true     | "api subdirectory file with depth 1"
        "docs/api/ref.md"   | "docs/api"  | 0     | true     | "immediate child of nested basePath"
        "docs/api/v1/spec.md" | "docs"    | 2     | true     | "deeply nested with sufficient depth"
        "docs/api/v1/spec.md" | "docs"    | 1     | false    | "deeply nested without sufficient depth"
    }

    def 'shouldIncludeAtDepth should handle realistic directory structure'() {
        given:
        def entries = [
            "/main.nf",
            "/nextflow.config",
            "/README.md",
            "/docs/guide.md",
            "/docs/api/index.md",
            "/docs/api/reference.md",
            "/src/process.nf",
            "/src/utils/helper.nf",
            "/test/test-data.csv"
        ]

        when: "listing root with depth 0"
        def rootDepth0 = entries.findAll { RepositoryProvider.shouldIncludeAtDepth(it, "/", 0) }

        then:
        rootDepth0.size() == 3
        rootDepth0.containsAll(["/main.nf", "/nextflow.config", "/README.md"])

        when: "listing root with depth 1"
        def rootDepth1 = entries.findAll { RepositoryProvider.shouldIncludeAtDepth(it, "/", 1) }

        then:
        rootDepth1.size() == 6
        rootDepth1.containsAll(["/main.nf", "/nextflow.config", "/README.md", "/docs/guide.md", "/src/process.nf", "/test/test-data.csv"])

        when: "listing docs with depth 0"
        def docsDepth0 = entries.findAll { RepositoryProvider.shouldIncludeAtDepth(it, "/docs", 0) }

        then:
        docsDepth0.size() == 1
        docsDepth0.contains("/docs/guide.md")

        when: "listing docs with depth 1"
        def docsDepth1 = entries.findAll { RepositoryProvider.shouldIncludeAtDepth(it, "/docs", 1) }

        then:
        docsDepth1.size() == 3
        docsDepth1.containsAll(["/docs/guide.md", "/docs/api/index.md", "/docs/api/reference.md"])
    }

}
