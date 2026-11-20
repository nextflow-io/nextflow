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
 * Offline tests for {@link BitbucketServerRepositoryProvider} against a WireMock fake of the Bitbucket Data Center REST API 1.0
 */
@Timeout(30)
class BitbucketServerRepositoryProviderWireMockTest extends ScmWireMockSpec {

    static final String REPO_PATH = '/rest/api/1.0/projects/acme/repos/hello'

    @Override
    protected String expectedAuth() { 'Bearer bbs-tkn' }

    private BitbucketServerRepositoryProvider provider(Map opts=[token: 'bbs-tkn']) {
        final config = new ProviderConfig('bbs', [server: url, platform: 'bitbucketserver'] + opts)
        return new BitbucketServerRepositoryProvider('scm/acme/hello', config)
    }

    def 'should return clone and repository url' () {
        given:
        stubJson(REPO_PATH, """
            {"slug": "hello", "name": "hello", "scmId": "git", "state": "AVAILABLE",
             "project": {"key": "ACME"},
             "links": {
                "clone": [
                    {"href": "ssh://git@localhost:7999/acme/hello.git", "name": "ssh"},
                    {"href": "$url/scm/acme/hello.git", "name": "http"}
                ],
                "self": [{"href": "$url/projects/ACME/repos/hello/browse"}]
             }}
            """)
        def repo = provider()

        expect:
        repo.getEndpointUrl() == "$url$REPO_PATH"
        repo.getCloneUrl() == "$url/scm/acme/hello.git"
        repo.getRepositoryUrl() == "$url/scm/acme/hello"
    }

    def 'should reject a non git repository' () {
        given:
        stubJson(REPO_PATH, '{"slug": "hello", "scmId": "hg", "links": {}}')

        when:
        provider().getCloneUrl()

        then:
        thrown(AbortOperationException)
    }

    def 'should read file content' () {
        given:
        stub("$REPO_PATH/raw/main.nf", ok('println "Hello world!"\n').withHeader('Content-Type', 'text/plain'))

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read binary content at revision with basic auth' () {
        given:
        stub("$REPO_PATH/raw/img/logo.png?at=feature%2Fx", ok().withBody(BINARY), 'Basic ' + 'me:secret'.bytes.encodeBase64().toString())

        expect:
        provider([user: 'me', password: 'secret']).setRevision('feature/x').readBytes('img/logo.png') == BINARY
    }

    def 'should list tags across pages' () {
        given:
        stubJson("$REPO_PATH/tags?start=0&limit=100", '''
            {"size": 1, "limit": 1, "isLastPage": false, "start": 0, "nextPageStart": 1, "values": [
                {"id": "refs/tags/v1.0", "displayId": "v1.0", "type": "TAG", "latestCommit": "t10", "hash": "h10"}
            ]}
            ''')
        stubJson("$REPO_PATH/tags?start=1&limit=100", '''
            {"size": 1, "limit": 1, "isLastPage": true, "start": 1, "values": [
                {"id": "refs/tags/v1.1", "displayId": "v1.1", "type": "TAG", "latestCommit": "t11", "hash": "h11"}
            ]}
            ''')

        expect:
        provider().getTags() == [
            new RepositoryProvider.TagInfo('v1.0', 't10'),
            new RepositoryProvider.TagInfo('v1.1', 't11') ]
    }

    def 'should list branches' () {
        given:
        stubJson("$REPO_PATH/branches?start=0&limit=100", '''
            {"size": 2, "limit": 100, "isLastPage": true, "start": 0, "values": [
                {"id": "refs/heads/master", "displayId": "master", "type": "BRANCH", "latestCommit": "aaa", "isDefault": true},
                {"id": "refs/heads/dev", "displayId": "dev", "type": "BRANCH", "latestCommit": "bbb", "isDefault": false}
            ]}
            ''')

        when:
        def result = provider().getBranches()

        then:
        result == [
            new RepositoryProvider.BranchInfo('master', 'aaa'),
            new RepositoryProvider.BranchInfo('dev', 'bbb') ]
        wireMock.verify(1, getRequestedFor(urlPathEqualTo("$REPO_PATH/branches")))
    }

    def 'should report missing repository' () {
        given:
        stubStatus("$REPO_PATH/raw/main.nf", 404, '{"errors": [{"message": "Repository acme/hello does not exist."}]}')

        when:
        provider().readText('main.nf')

        then:
        thrown(AbortOperationException)
    }
}
