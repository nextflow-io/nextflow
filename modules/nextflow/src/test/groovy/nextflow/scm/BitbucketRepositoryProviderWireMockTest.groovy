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
import spock.lang.Timeout

/**
 * Offline tests for {@link BitbucketRepositoryProvider} against a WireMock fake of the Bitbucket Cloud REST API 2.0
 */
@Timeout(30)
class BitbucketRepositoryProviderWireMockTest extends ScmWireMockSpec {

    static final String REPO_PATH = '/2.0/repositories/acme/hello'

    @Override
    protected String expectedAuth() {
        'Basic ' + 'me:app-pass'.bytes.encodeBase64().toString()
    }

    private BitbucketRepositoryProvider provider(Map opts=[user: 'me', password: 'app-pass']) {
        final config = new ProviderConfig('bitbucket', [server: url, endpoint: url] + opts)
        return new BitbucketRepositoryProvider('acme/hello', config)
    }

    private void stubRepo(String scm='git') {
        stubJson(REPO_PATH, """
            {"type": "repository", "full_name": "acme/hello", "scm": "$scm",
             "mainbranch": {"type": "branch", "name": "main"},
             "links": {"clone": [
                {"name": "https", "href": "https://me@bitbucket.org/acme/hello.git"},
                {"name": "ssh", "href": "git@bitbucket.org:acme/hello.git"}
             ]}}
            """)
    }

    def 'should return clone url, main branch and repository url' () {
        given:
        stubRepo()
        def repo = provider()

        expect:
        repo.getEndpointUrl() == "$url$REPO_PATH"
        repo.getCloneUrl() == 'https://me@bitbucket.org/acme/hello.git'
        repo.getMainBranch() == 'main'
        repo.getRepositoryUrl() == "$url/acme/hello"
    }

    def 'should reject a non git repository' () {
        given:
        stubRepo('hg')

        when:
        provider().getCloneUrl()

        then:
        thrown(AbortOperationException)
    }

    def 'should read file content from the main branch' () {
        given:
        stubRepo()
        stub("$REPO_PATH/src/main/main.nf", ok('println "Hello world!"\n').withHeader('Content-Type', 'text/plain'))

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read binary content at revision' () {
        given:
        stub("$REPO_PATH/src/v1.0/img/logo.png", ok().withBody(BINARY))

        when:
        def result = provider().setRevision('v1.0').readBytes('img/logo.png')

        then:
        result == BINARY
        // a revision without slashes is used as is
        wireMock.verify(0, getRequestedFor(urlPathMatching("$REPO_PATH/refs/.*")))
    }

    def 'should resolve a #kind name with a slash to its commit' () {
        given:
        if( kind == 'tag' )
            stubStatus("$REPO_PATH/refs/branches/$revision", 404, '{"type": "error"}')
        stubJson("$REPO_PATH/refs/$refs/$revision", """
            {"type": "$kind", "name": "$revision", "target": {"type": "commit", "hash": "$commit"}}
            """)
        stub("$REPO_PATH/src/$commit/main.nf", ok(kind))

        expect:
        provider().setRevision(revision).readText('main.nf') == kind

        where:
        kind     | refs       | revision      | commit
        'branch' | 'branches' | 'feature/x'   | 'abc123'
        'tag'    | 'tags'     | 'release/1.0' | 'def456'
    }

    def 'should authenticate with api token' () {
        given:
        stub("$REPO_PATH/src/dev/main.nf", ok('hello'), 'Basic ' + 'me@acme.com:api-tkn'.bytes.encodeBase64().toString())

        expect:
        provider([user: 'me@acme.com', token: 'api-tkn']).setRevision('dev').readText('main.nf') == 'hello'
    }

    def 'should list tags following the next page link' () {
        given:
        stubJson("$REPO_PATH/refs/tags", """
            {"pagelen": 1, "page": 1, "size": 2, "next": "$url$REPO_PATH/refs/tags?page=2",
             "values": [{"type": "tag", "name": "v1.0", "target": {"type": "commit", "hash": "t10"}}]}
            """)
        stubJson("$REPO_PATH/refs/tags?page=2", '''
            {"pagelen": 1, "page": 2, "size": 2,
             "values": [{"type": "tag", "name": "v1.1", "target": {"type": "commit", "hash": "t11"}}]}
            ''')

        when:
        def result = provider().getTags()

        then:
        result == [
            new RepositoryProvider.TagInfo('v1.0', 't10'),
            new RepositoryProvider.TagInfo('v1.1', 't11') ]
        wireMock.verify(2, getRequestedFor(urlPathEqualTo("$REPO_PATH/refs/tags")))
    }

    def 'should list branches' () {
        given:
        stubJson("$REPO_PATH/refs/branches", '''
            {"pagelen": 10, "page": 1, "size": 2, "values": [
                {"type": "branch", "name": "main", "target": {"type": "commit", "hash": "aaa"}},
                {"type": "branch", "name": "dev", "target": {"type": "commit", "hash": "bbb"}}
            ]}
            ''')

        expect:
        provider().getBranches() == [
            new RepositoryProvider.BranchInfo('main', 'aaa'),
            new RepositoryProvider.BranchInfo('dev', 'bbb') ]
    }

    def 'should list root directory' () {
        given:
        stubRepo()
        stubJson("$REPO_PATH/src/main/", '''
            {"pagelen": 10, "page": 1, "values": [
                {"type": "commit_file", "path": "main.nf", "size": 42, "commit": {"type": "commit", "hash": "c1"}},
                {"type": "commit_directory", "path": "docs", "commit": {"type": "commit", "hash": "c1"}}
            ]}
            ''')

        when:
        def result = provider().listDirectory('/', 1)

        then:
        result == [
            new RepositoryProvider.RepositoryEntry('docs', '/docs', RepositoryProvider.EntryType.DIRECTORY, 'c1', null),
            new RepositoryProvider.RepositoryEntry('main.nf', '/main.nf', RepositoryProvider.EntryType.FILE, 'c1', 42) ]
    }

    def 'should list subdirectory at revision' () {
        given:
        stubJson("$REPO_PATH/src/dev/docs", '''
            {"values": [
                {"type": "commit_file", "path": "docs/index.md", "size": 10, "commit": {"hash": "c2"}},
                {"type": "commit_directory", "path": "docs/guide", "commit": {"hash": "c2"}}
            ]}
            ''')

        when:
        def result = provider().setRevision('dev').listDirectory('docs', 1)

        then:
        result*.path == ['/docs/guide', '/docs/index.md']
        result*.type == [RepositoryProvider.EntryType.DIRECTORY, RepositoryProvider.EntryType.FILE]
    }

    def 'should report directory listing errors as unsupported' () {
        given:
        stubStatus("$REPO_PATH/src/dev/missing", 404, '{"type": "error"}')

        when:
        provider().setRevision('dev').listDirectory('missing', 1)

        then:
        def e = thrown(UnsupportedOperationException)
        e.cause instanceof AbortOperationException
    }
}
