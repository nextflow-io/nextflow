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

import nextflow.exception.AbortOperationException
import spock.lang.PendingFeature
import spock.lang.Timeout

/**
 * Offline tests for {@link GiteaRepositoryProvider} against a WireMock fake of the Gitea REST API v1
 */
@Timeout(30)
class GiteaRepositoryProviderWireMockTest extends ScmWireMockSpec {

    static final String REPO_PATH = '/api/v1/repos/acme/hello'

    @Override
    protected String expectedAuth() { 'token gt-tkn' }

    private GiteaRepositoryProvider provider(Map opts=[token: 'gt-tkn']) {
        // the endpoint defaults to `<server>/api/v1`
        final config = new ProviderConfig('gitea', [server: url] + opts)
        return new GiteaRepositoryProvider('acme/hello', config)
    }

    def 'should return clone and repository url' () {
        given:
        stubJson(REPO_PATH, """
            {"id": 3, "name": "hello", "full_name": "acme/hello", "default_branch": "main",
             "clone_url": "$url/acme/hello.git", "ssh_url": "git@localhost:acme/hello.git",
             "html_url": "$url/acme/hello"}
            """)
        def repo = provider()

        expect:
        repo.getEndpointUrl() == "$url$REPO_PATH"
        repo.getCloneUrl() == "$url/acme/hello.git"
        repo.getRepositoryUrl() == "$url/acme/hello"
    }

    def 'should read file content' () {
        given:
        stub("$REPO_PATH/raw/main.nf", ok('println "Hello world!"\n').withHeader('Content-Type', 'text/plain'))

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read binary content at revision' () {
        given:
        stub("$REPO_PATH/raw/img/logo.png?ref=feature%2Fx", ok().withBody(BINARY))

        expect:
        provider().setRevision('feature/x').readBytes('img/logo.png') == BINARY
    }

    def 'should send no auth header with user and password' () {
        given:
        // known gap: Gitea accepts basic auth, but only the token is sent
        stub("$REPO_PATH/raw/main.nf", ok('hello'), NO_AUTH)

        expect:
        provider([user: 'me', password: 'secret']).readText('main.nf') == 'hello'
    }

    @PendingFeature(reason = 'user and password are not sent as basic auth')
    def 'should send basic auth with user and password' () {
        given:
        stub("$REPO_PATH/raw/main.nf", ok('hello'), 'Basic ' + 'me:secret'.bytes.encodeBase64().toString())

        expect:
        provider([user: 'me', password: 'secret']).readText('main.nf') == 'hello'
    }

    def 'should list branches and tags' () {
        given:
        stubJson("$REPO_PATH/branches?page=1", '''
            [{"name": "main", "commit": {"id": "aaa", "message": "x"}, "protected": false},
             {"name": "dev", "commit": {"id": "bbb", "message": "y"}, "protected": false}]
            ''')
        stubJson("$REPO_PATH/branches?page=2", '[]')
        stubJson("$REPO_PATH/tags?page=1", '''
            [{"name": "v1.0", "id": "tag-obj", "commit": {"sha": "t10", "url": "x"}}]
            ''')
        stubJson("$REPO_PATH/tags?page=2", '[]')
        def repo = provider()

        expect:
        repo.getBranches() == [
            new RepositoryProvider.BranchInfo('main', 'aaa'),
            new RepositoryProvider.BranchInfo('dev', 'bbb') ]
        repo.getTags() == [ new RepositoryProvider.TagInfo('v1.0', 't10') ]
    }

    def 'should list root directory' () {
        given:
        // known quirk: the listing defaults to master instead of the repo default branch
        stubJson("$REPO_PATH/contents?ref=master", '''
            [{"name": "main.nf", "path": "main.nf", "sha": "f1", "type": "file", "size": 42},
             {"name": "docs", "path": "docs", "sha": "d1", "type": "dir", "size": 0}]
            ''')

        when:
        def result = provider().listDirectory('/', 1)

        then:
        result == [
            new RepositoryProvider.RepositoryEntry('docs', '/docs', RepositoryProvider.EntryType.DIRECTORY, 'd1', 0),
            new RepositoryProvider.RepositoryEntry('main.nf', '/main.nf', RepositoryProvider.EntryType.FILE, 'f1', 42) ]
    }

    def 'should list subdirectory recursively at revision' () {
        given:
        stubJson("$REPO_PATH/contents/docs?ref=dev", '''
            [{"name": "index.md", "path": "docs/index.md", "sha": "f2", "type": "file", "size": 10},
             {"name": "guide", "path": "docs/guide", "sha": "d2", "type": "dir", "size": 0}]
            ''')
        stubJson("$REPO_PATH/contents/docs/guide?ref=dev", '''
            [{"name": "intro.md", "path": "docs/guide/intro.md", "sha": "f3", "type": "file", "size": 20}]
            ''')

        when:
        def result = provider().setRevision('dev').listDirectory('docs', 2)

        then:
        result*.path == ['/docs/guide', '/docs/index.md', '/docs/guide/intro.md']
    }

    def 'should report directory listing errors as unsupported' () {
        given:
        stubStatus("$REPO_PATH/contents/missing?ref=master", 404, '{"message": "GetContentsOrList"}')

        when:
        provider().listDirectory('missing', 1)

        then:
        def e = thrown(UnsupportedOperationException)
        e.cause instanceof AbortOperationException
    }
}
