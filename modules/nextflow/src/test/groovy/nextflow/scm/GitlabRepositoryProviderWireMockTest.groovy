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
 * Offline tests for {@link GitlabRepositoryProvider} against a WireMock fake of the GitLab REST API v4
 */
@Timeout(30)
class GitlabRepositoryProviderWireMockTest extends ScmWireMockSpec {

    static final String PROJECT_PATH = '/api/v4/projects/acme%2Fhello'

    @Override
    protected String authHeader() { 'PRIVATE-TOKEN' }

    @Override
    protected String expectedAuth() { 'glpat-123' }

    private GitlabRepositoryProvider provider(Map opts=[token: 'glpat-123']) {
        final config = new ProviderConfig('gitlab', [server: url] + opts)
        return new GitlabRepositoryProvider('acme/hello', config)
    }

    private void stubProject() {
        stubJson(PROJECT_PATH, """
            {"id": 7, "name": "hello", "path_with_namespace": "acme/hello", "default_branch": "main",
             "http_url_to_repo": "$url/acme/hello.git", "ssh_url_to_repo": "git@localhost:acme/hello.git",
             "web_url": "$url/acme/hello"}
            """)
    }

    def 'should return clone url, default branch and repository url' () {
        given:
        stubProject()
        def repo = provider()

        expect:
        repo.getEndpointUrl() == "$url$PROJECT_PATH"
        repo.getCloneUrl() == "$url/acme/hello.git"
        repo.getDefaultBranch() == 'main'
        repo.getRepositoryUrl() == "$url/acme/hello"
    }

    def 'should fall back to master when the default branch is not available' () {
        given:
        stubJson(PROJECT_PATH, '{"id": 7, "name": "hello"}')

        expect:
        provider().getDefaultBranch() == 'master'
    }

    def 'should read file content from the default branch' () {
        given:
        stubProject()
        stubJson("$PROJECT_PATH/repository/files/main.nf?ref=main", """
            {"file_name": "main.nf", "file_path": "main.nf", "size": 23, "encoding": "base64",
             "ref": "main", "blob_id": "b1", "commit_id": "c1",
             "content": "${'println "Hello world!"\n'.bytes.encodeBase64()}"}
            """)

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read binary content at revision with encoded path' () {
        given:
        stubJson("$PROJECT_PATH/repository/files/img%2Flogo.png?ref=feature%2Fx", """
            {"file_path": "img/logo.png", "encoding": "base64", "content": "${BINARY.encodeBase64()}"}
            """)

        when:
        def result = provider().setRevision('feature/x').readBytes('/img/logo.png')

        then:
        result == BINARY
        // the revision is given, so no default branch lookup
        wireMock.verify(0, getRequestedFor(urlEqualTo(PROJECT_PATH)))
    }

    def 'should send the password as private token' () {
        given:
        stubJson("$PROJECT_PATH/repository/files/main.nf?ref=dev", '{"content": "aGVsbG8="}', 'secret')

        expect:
        provider([user: 'me', password: 'secret']).setRevision('dev').readText('main.nf') == 'hello'
    }

    def 'should list branches and tags' () {
        given:
        stubJson("$PROJECT_PATH/repository/branches?page=1", '''
            [{"name": "main", "default": true, "commit": {"id": "aaa", "short_id": "aaa"}},
             {"name": "dev", "default": false, "commit": {"id": "bbb", "short_id": "bbb"}}]
            ''')
        stubJson("$PROJECT_PATH/repository/branches?page=2", '[]')
        stubJson("$PROJECT_PATH/repository/tags?page=1", '''
            [{"name": "v1.0", "message": "", "target": "t10", "commit": {"id": "c10"}}]
            ''')
        stubJson("$PROJECT_PATH/repository/tags?page=2", '[]')
        def repo = provider()

        expect:
        repo.getBranches() == [
            new RepositoryProvider.BranchInfo('main', 'aaa'),
            new RepositoryProvider.BranchInfo('dev', 'bbb') ]
        repo.getTags() == [ new RepositoryProvider.TagInfo('v1.0', 'c10') ]
    }

    def 'should list root directory' () {
        given:
        stubJson("$PROJECT_PATH/repository/tree?ref=dev", '''
            [{"id": "f1", "name": "main.nf", "type": "blob", "path": "main.nf", "mode": "100644"},
             {"id": "d1", "name": "docs", "type": "tree", "path": "docs", "mode": "040000"}]
            ''')

        when:
        def result = provider().setRevision('dev').listDirectory('/', 1)

        then:
        result == [
            new RepositoryProvider.RepositoryEntry('docs', '/docs', RepositoryProvider.EntryType.DIRECTORY, 'd1', null),
            new RepositoryProvider.RepositoryEntry('main.nf', '/main.nf', RepositoryProvider.EntryType.FILE, 'f1', null) ]
    }

    def 'should list subdirectory recursively' () {
        given:
        stubProject()
        stubJson("$PROJECT_PATH/repository/tree?ref=main&path=docs%2Fguide&recursive=true", '''
            [{"id": "f2", "name": "intro.md", "type": "blob", "path": "docs/guide/intro.md"},
             {"id": "d2", "name": "deep", "type": "tree", "path": "docs/guide/deep"},
             {"id": "f3", "name": "more.md", "type": "blob", "path": "docs/guide/deep/more.md"},
             {"id": "f4", "name": "x.md", "type": "blob", "path": "docs/guide/deep/deeper/x.md"},
             {"id": "f5", "name": "y.md", "type": "blob", "path": "docs/guide/deep/deeper/deepest/y.md"}]
            ''')

        when:
        def result = provider().listDirectory('docs/guide', 2)

        then:
        result*.path == ['/docs/guide/deep', '/docs/guide/intro.md', '/docs/guide/deep/more.md', '/docs/guide/deep/deeper/x.md']
    }

    def 'should report missing project' () {
        given:
        stubStatus(PROJECT_PATH, 404, '{"message": "404 Project Not Found"}')

        when:
        provider().getCloneUrl()

        then:
        thrown(AbortOperationException)
    }
}
