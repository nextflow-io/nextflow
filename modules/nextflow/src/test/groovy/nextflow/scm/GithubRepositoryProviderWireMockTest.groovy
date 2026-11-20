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

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.exception.RateLimitExceededException
import spock.lang.Timeout

/**
 * Offline tests for {@link GithubRepositoryProvider} against a WireMock fake of the GitHub REST API
 */
@Timeout(30)
class GithubRepositoryProviderWireMockTest extends ScmWireMockSpec {

    @Override
    protected String expectedAuth() {
        'Basic ' + 'tkn123:x-oauth-basic'.bytes.encodeBase64().toString()
    }

    private GithubRepositoryProvider provider(Map opts=[token: 'tkn123']) {
        final config = new ProviderConfig('github', [server: url, endpoint: url] + opts)
        return new GithubRepositoryProvider('acme/hello', config)
    }

    def 'should read file content' () {
        given:
        // the contents API wraps base64 content every 60 chars
        final content = 'println "Hello world!"\n'.bytes.encodeBase64().toString()
        stubJson('/repos/acme/hello/contents/main.nf', """
            {"type": "file", "encoding": "base64", "name": "main.nf", "path": "main.nf", "sha": "abc",
             "content": "${content.substring(0,10)}\\n${content.substring(10)}\\n"}
            """)

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read binary content at revision' () {
        given:
        stubJson('/repos/acme/hello/contents/img/logo.png?ref=feature%2Fx', """
            {"type": "file", "encoding": "base64", "content": "${BINARY.encodeBase64()}"}
            """)

        expect:
        provider().setRevision('feature/x').readBytes('img/logo.png') == BINARY
    }

    def 'should not send auth header without credentials' () {
        given:
        stubJson('/repos/acme/hello/contents/main.nf', '{"content": "aGVsbG8="}', NO_AUTH)

        expect:
        provider([:]).readText('main.nf') == 'hello'
    }

    def 'should authenticate with the GITHUB_TOKEN env variable' () {
        given:
        SysEnv.push([GITHUB_TOKEN: 'envtkn'])
        stubJson('/repos/acme/hello/contents/main.nf', '{"content": "aGVsbG8="}',
            'Basic ' + 'envtkn:x-oauth-basic'.bytes.encodeBase64().toString())

        expect:
        provider([:]).readText('main.nf') == 'hello'

        cleanup:
        SysEnv.pop()
    }

    def 'should return clone and repository url' () {
        given:
        stubJson('/repos/acme/hello', """
            {"id": 1, "name": "hello", "full_name": "acme/hello", "default_branch": "main",
             "clone_url": "https://github.com/acme/hello.git", "html_url": "https://github.com/acme/hello"}
            """)
        def repo = provider()

        expect:
        repo.getCloneUrl() == 'https://github.com/acme/hello.git'
        repo.getRepositoryUrl() == "$url/acme/hello"
        repo.getEndpointUrl() == "$url/repos/acme/hello"
    }

    def 'should list branches across pages' () {
        given:
        stubJson('/repos/acme/hello/branches?page=1', '''
            [{"name": "main", "commit": {"sha": "aaa", "url": "x"}, "protected": false},
             {"name": "dev", "commit": {"sha": "bbb", "url": "x"}, "protected": false}]
            ''')
        stubJson('/repos/acme/hello/branches?page=2', '''
            [{"name": "feature/x", "commit": {"sha": "ccc", "url": "x"}, "protected": false}]
            ''')
        stubJson('/repos/acme/hello/branches?page=3', '[]')

        when:
        def result = provider().getBranches()

        then:
        result == [
            new RepositoryProvider.BranchInfo('main', 'aaa'),
            new RepositoryProvider.BranchInfo('dev', 'bbb'),
            new RepositoryProvider.BranchInfo('feature/x', 'ccc') ]
        wireMock.verify(3, getRequestedFor(urlPathEqualTo('/repos/acme/hello/branches')))
    }

    def 'should list tags' () {
        given:
        stubJson('/repos/acme/hello/tags?page=1', '''
            [{"name": "v1.1", "commit": {"sha": "t11", "url": "x"}, "zipball_url": "x", "tarball_url": "x"},
             {"name": "v1.0", "commit": {"sha": "t10", "url": "x"}, "zipball_url": "x", "tarball_url": "x"}]
            ''')
        stubJson('/repos/acme/hello/tags?page=2', '[]')

        expect:
        provider().getTags() == [
            new RepositoryProvider.TagInfo('v1.1', 't11'),
            new RepositoryProvider.TagInfo('v1.0', 't10') ]
    }

    def 'should list root directory falling back to the master branch' () {
        given:
        // known quirk: the main-then-master fallback ignores the repo default_branch
        stubStatus('/repos/acme/hello/git/refs/heads/main', 404, '{"message": "Not Found"}')
        stubJson('/repos/acme/hello/git/refs/heads/master', '''
            {"ref": "refs/heads/master", "object": {"sha": "c0ffee", "type": "commit"}}
            ''')
        stubJson('/repos/acme/hello/git/commits/c0ffee', '''
            {"sha": "c0ffee", "tree": {"sha": "root-tree", "url": "x"}, "message": "init"}
            ''')
        stubJson('/repos/acme/hello/git/trees/root-tree', '''
            {"sha": "root-tree", "truncated": false, "tree": [
                {"path": "main.nf", "mode": "100644", "type": "blob", "sha": "f1", "size": 42},
                {"path": "docs", "mode": "040000", "type": "tree", "sha": "docs-tree"}
            ]}
            ''')

        when:
        def result = provider().listDirectory('/', 1)

        then:
        result == [
            new RepositoryProvider.RepositoryEntry('docs', '/docs', RepositoryProvider.EntryType.DIRECTORY, 'docs-tree', null),
            new RepositoryProvider.RepositoryEntry('main.nf', '/main.nf', RepositoryProvider.EntryType.FILE, 'f1', 42) ]
    }

    def 'should list subdirectory recursively at revision' () {
        given:
        stubJson('/repos/acme/hello/git/refs/heads/dev', '''
            {"ref": "refs/heads/dev", "object": {"sha": "d3v", "type": "commit"}}
            ''')
        stubJson('/repos/acme/hello/git/commits/d3v', '''
            {"sha": "d3v", "tree": {"sha": "root-tree"}}
            ''')
        stubJson('/repos/acme/hello/git/trees/root-tree', '''
            {"sha": "root-tree", "tree": [
                {"path": "main.nf", "type": "blob", "sha": "f1", "size": 42},
                {"path": "docs", "type": "tree", "sha": "docs-tree"}
            ]}
            ''')
        // a recursive tree includes the intermediate tree entries
        stubJson('/repos/acme/hello/git/trees/docs-tree?recursive=1', '''
            {"sha": "docs-tree", "truncated": false, "tree": [
                {"path": "index.md", "type": "blob", "sha": "f2", "size": 10},
                {"path": "guide", "type": "tree", "sha": "guide-tree"},
                {"path": "guide/intro.md", "type": "blob", "sha": "f3", "size": 20},
                {"path": "guide/deep", "type": "tree", "sha": "deep-tree"},
                {"path": "guide/deep/more.md", "type": "blob", "sha": "f4", "size": 30}
            ]}
            ''')

        when:
        def result = provider().setRevision('dev').listDirectory('docs', 2)

        then:
        result*.path == ['/docs/guide/deep', '/docs/guide', '/docs/index.md', '/docs/guide/intro.md', '/docs/guide/deep/more.md']
        result.findAll { it.type == RepositoryProvider.EntryType.DIRECTORY }*.path == ['/docs/guide/deep', '/docs/guide']
    }

    def 'should report error status #status' () {
        given:
        final response = aResponse().withStatus(status).withBody(body)
        if( rateLimited )
            response.withHeader('X-RateLimit-Remaining', '0')
        stub('/repos/acme/hello/contents/main.nf', response, credentials ? expectedAuth() : NO_AUTH)

        when:
        provider(credentials ? [token: 'tkn123'] : [:]).readText('main.nf')

        then:
        def e = thrown(Exception)
        e.class == exception
        e.message.startsWith(message)

        where:
        status | rateLimited | credentials | body                                     | exception                  | message
        404    | false       | true        | '{"message": "Not Found"}'               | AbortOperationException    | 'Remote resource not found: '
        403    | true        | false       | '{"message": "API rate limit exceeded"}' | RateLimitExceededException | 'API rate limit exceeded'
        401    | false       | true        | '{"message": "Bad credentials"}'         | AbortOperationException    | 'Not authorized'
    }
}
