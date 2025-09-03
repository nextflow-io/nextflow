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
 */

package nextflow.scm

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.channels.UnresolvedAddressException
import javax.net.ssl.SSLSession

import nextflow.SysEnv
import nextflow.exception.HttpResponseLengthExceedException
import nextflow.util.RetryConfig
import spock.lang.Specification
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

}
