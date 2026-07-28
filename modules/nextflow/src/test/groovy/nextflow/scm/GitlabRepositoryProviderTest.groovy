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
import java.net.http.HttpHeaders
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import javax.net.ssl.SSLSession

import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@IgnoreIf({System.getenv('NXF_SMOKE')})
class GitlabRepositoryProviderTest extends Specification {

    def 'should return repo url' () {
        expect:
        new GitlabRepositoryProvider('pditommaso/hello').getEndpointUrl() == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello'
    }

    def 'should return project URL' () {
        expect:
        new GitlabRepositoryProvider('pditommaso/hello').getRepositoryUrl() == 'https://gitlab.com/pditommaso/hello'
    }

    @Unroll
    def 'should validate hasCredentials' () {
        given:
        def provider = new GitlabRepositoryProvider('pditommaso/tutorial', CONFIG)

        expect:
        provider.hasCredentials() == EXPECTED

        where:
        EXPECTED    | CONFIG
        false       | new ProviderConfig('gitlab')
        false       | new ProviderConfig('gitlab').setUser('foo')
        true        | new ProviderConfig('gitlab').setUser('foo').setPassword('bar')
        true        | new ProviderConfig('gitlab').setToken('xyz')
    }

    @Unroll
    def 'should return git credentials' () {
        given:
        def provider = new GitlabRepositoryProvider('pditommaso/tutorial', CONFIG)

        when:
        def credentials = provider.getGitCredentials()

        then:
        credentials != null

        where:
        CONFIG                                                                    | _
        new ProviderConfig('gitlab').setUser('foo').setPassword('bar')            | _
        new ProviderConfig('gitlab').setUser('foo').setToken('xyz')               | _
        new ProviderConfig('gitlab').setUser('foo').setPassword('bar').setToken('xyz') | _
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should return clone url'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)

        when:
        def url = new GitlabRepositoryProvider('pditommaso/hello', config).getCloneUrl()
        then:
        url == 'https://gitlab.com/pditommaso/hello.git'
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should read file content'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)

        when:
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)
        def result = repo.readText('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

        when:
        repo = new GitlabRepositoryProvider('pditommaso/hello', config)
        repo.setRevision('test/branch+with&special-chars')
        result = repo.readText('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should read binary content'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)
        and:
        def DATA = this.class.getResourceAsStream('/test-asset.bin').bytes

        when:
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)
        def result = repo.readBytes('test/test-asset.bin')
        then:
        result == DATA
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should return default branch' () {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)

        when:
        def provider = new GitlabRepositoryProvider('pditommaso/hello', config)
        then:
        provider.getDefaultBranch() == 'master'
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should return content URL' () {
        given:
        def (user, pwd, token) = System.getenv('NXF_GITLAB_ACCESS_TOKEN').tokenize(':')
        String CONFIG = """
        providers {
            mygitlab {
                server = 'https://gitlab.com'
                endpoint = 'https://gitlab.com'
                platform = 'gitlab'
                user = '$user'
                password = '$token' // NOTE: Gitlab token can be used in place of the password
            }
        }
        """

        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('github', config.providers.mygitlab as ConfigObject)

        expect:
        new GitlabRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/main.nf?ref=master'

        and:
        new GitlabRepositoryProvider('pditommaso/hello', obj)
                .setRevision('the-commit-id')
                .getContentUrl('main.nf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/main.nf?ref=the-commit-id'

        and:
        new GitlabRepositoryProvider('pditommaso/hello', obj)
                .setRevision('test/branch+with&strangecharacters')
                .getContentUrl('main.nf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/main.nf?ref=test%2Fbranch%2Bwith%26strangecharacters'

        and:
        new GitlabRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('conf/extra.conf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/conf%2Fextra.conf?ref=master'


        and: // should strip leading slashes
        new GitlabRepositoryProvider('pditommaso/hello', obj)
            .getContentUrl('/main.nf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/main.nf?ref=master'

        and:
        new GitlabRepositoryProvider('pditommaso/hello', obj)
            .getContentUrl('//conf/extra.conf') == 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/files/conf%2Fextra.conf?ref=master'

    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should list root directory contents'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)

        when:
        def entries = repo.listDirectory("/", 1)

        then:
        entries.size() > 0
        and:
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.any { it.name == 'test' && it.type == RepositoryProvider.EntryType.DIRECTORY }
        and:
        // Should NOT include nested files for depth=1
        !entries.any { it.path == '/test/test-asset.bin' }
        and:
        entries.every { it.path && it.sha }
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should list subdirectory contents'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)

        when:
        def entries = repo.listDirectory("/test", 1)

        then:
        entries.size() > 0
        entries.any { it.name == 'test-asset.bin' && it.path=='/test/test-asset.bin' && it.type == RepositoryProvider.EntryType.FILE }
        entries.every { it.path.startsWith('/test/') }
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should list directory contents recursively'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)

        when:
        def entries = repo.listDirectory("/", 10)

        then:
        entries.size() > 0
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.any { it.name == 'test-asset.bin' && it.type == RepositoryProvider.EntryType.FILE }
        entries.every { it.path && it.sha }
    }

    @Requires({System.getenv('NXF_GITLAB_ACCESS_TOKEN')})
    def 'should list directory contents with depth 2'() {
        given:
        def token = System.getenv('NXF_GITLAB_ACCESS_TOKEN')
        def config = new ProviderConfig('gitlab').setAuth(token)
        def repo = new GitlabRepositoryProvider('pditommaso/hello', config)

        when:
        def entries = repo.listDirectory("/", 2)

        then:
        entries.size() > 0
        // Should include immediate children (depth 1)
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.any { it.name == 'test' && it.type == RepositoryProvider.EntryType.DIRECTORY }
        // Should include nested files (depth 2)
        entries.any { it.name == 'test-asset.bin' && it.path.contains('/test/') }
        entries.every { it.path && it.sha }
    }

    def 'should follow GitLab pagination links when listing branches' () {
        given:
        def provider = Spy(GitlabRepositoryProvider, constructorArgs: ['pditommaso/hello', new ProviderConfig('gitlab')])
        and:
        def first = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/branches?per_page=100'
        def second = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/branches?per_page=100&page=2'
        def last = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/branches?per_page=100&page=3'

        when:
        def branches = provider.getBranches()

        then:
        1 * provider.invokeResponse(first) >> response(
            first,
            '[{"name":"main","commit":{"id":"aaa"}}]',
            "<${last}>; rel=\"last\", <${second}>; type=\"application/json\"; rel=\"prev next\""
        )
        1 * provider.invokeResponse(second) >> response(
            second,
            '[{"name":"develop","commit":{"id":"bbb"}}]'
        )
        and:
        branches.name == ['main', 'develop']
    }

    def 'should request the max page size when listing tags' () {
        given:
        def provider = Spy(GitlabRepositoryProvider, constructorArgs: ['pditommaso/hello', new ProviderConfig('gitlab')])
        and:
        def url = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/tags?per_page=100'

        when:
        def tags = provider.getTags()

        then:
        1 * provider.invokeResponse(url) >> response(
            url,
            '[{"name":"v1.0","commit":{"id":"aaa"}}]'
        )
        and:
        tags.name == ['v1.0']
    }

    def 'should reject a cross-origin GitLab pagination link' () {
        given:
        def provider = Spy(GitlabRepositoryProvider, constructorArgs: ['pditommaso/hello', new ProviderConfig('gitlab')])
        and:
        def url = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/branches?per_page=100'

        when:
        provider.getBranches()

        then:
        1 * provider.invokeResponse(url) >> response(
            url,
            '[{"name":"main","commit":{"id":"aaa"}}]',
            '<https://example.com/api/v4/projects/1/repository/branches?page=2>; rel="next"'
        )
        def error = thrown(IOException)
        error.message == 'Invalid GitLab pagination URL: https://example.com/api/v4/projects/1/repository/branches?page=2'
    }

    def 'should reject a GitLab pagination link cycle' () {
        given:
        def provider = Spy(GitlabRepositoryProvider, constructorArgs: ['pditommaso/hello', new ProviderConfig('gitlab')])
        and:
        def url = 'https://gitlab.com/api/v4/projects/pditommaso%2Fhello/repository/branches?per_page=100'

        when:
        provider.getBranches()

        then:
        1 * provider.invokeResponse(url) >> response(
            url,
            '[{"name":"main","commit":{"id":"aaa"}}]',
            "<${url}>; rel=\"next\""
        )
        def error = thrown(IOException)
        error.message == "Invalid GitLab pagination link cycle detected: $url"
    }

    private static HttpResponse<byte[]> response(String url, String responseBody, String link=null) {
        return new HttpResponse<byte[]>() {
            @Override
            int statusCode() { 200 }

            @Override
            HttpRequest request() { null }

            @Override
            Optional<HttpResponse<byte[]>> previousResponse() { Optional.empty() }

            @Override
            HttpHeaders headers() {
                final values = link ? ['Link': [link]] : [:]
                HttpHeaders.of(values, (a, b) -> true)
            }

            @Override
            byte[] body() { responseBody.bytes }

            @Override
            Optional<SSLSession> sslSession() { Optional.empty() }

            @Override
            URI uri() { new URI(url) }

            @Override
            HttpClient.Version version() { HttpClient.Version.HTTP_1_1 }
        }
    }
}
