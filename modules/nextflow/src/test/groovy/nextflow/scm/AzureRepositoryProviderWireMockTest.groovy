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
 * Offline tests for {@link AzureRepositoryProvider} against a WireMock fake of the Azure DevOps Git REST API 6.0
 */
@Timeout(30)
class AzureRepositoryProviderWireMockTest extends ScmWireMockSpec {

    static final String REPO_PATH = '/acme/proj/_apis/git/repositories/hello'

    static final String COMMIT = 'a' * 40

    static final String TAG_OBJECT = 'b' * 40

    static final String TAG_COMMIT = 'c' * 40

    @Override
    protected String expectedAuth() {
        'Basic ' + 'me:pat-123'.bytes.encodeBase64().toString()
    }

    private AzureRepositoryProvider provider(Map opts=[user: 'me', password: 'pat-123']) {
        final config = new ProviderConfig('azurerepos', [server: url, endpoint: url] + opts)
        return new AzureRepositoryProvider('acme/proj/hello', config)
    }

    /*
     * v1.0 is lightweight and v1.1 is annotated: the objectId of an annotated tag is the
     * tag object, and the commit is returned as peeledObjectId only when peelTags=true is sent
     */
    private void stubTags() {
        final tags = get(urlPathEqualTo("$REPO_PATH/refs")).withQueryParam('filter', equalTo('tags'))
        wireMock.stubFor(withAuth(tags.withQueryParam('peelTags', absent())).willReturn(okJson("""
            {"count": 2, "value": [
                {"name": "refs/tags/v1.0", "objectId": "t10", "creator": {"displayName": "x"}},
                {"name": "refs/tags/v1.1", "objectId": "$TAG_OBJECT", "creator": {"displayName": "x"}}
            ]}
            """)))
        final peeled = get(urlPathEqualTo("$REPO_PATH/refs")).withQueryParam('filter', equalTo('tags'))
        wireMock.stubFor(withAuth(peeled.withQueryParam('peelTags', equalTo('true'))).willReturn(okJson("""
            {"count": 2, "value": [
                {"name": "refs/tags/v1.0", "objectId": "t10", "creator": {"displayName": "x"}},
                {"name": "refs/tags/v1.1", "objectId": "$TAG_OBJECT", "peeledObjectId": "$TAG_COMMIT", "creator": {"displayName": "x"}}
            ]}
            """)))
    }

    def 'should return endpoint, clone and repository url' () {
        given:
        def repo = provider()

        expect:
        repo.getEndpointUrl() == "$url$REPO_PATH"
        repo.getRepositoryUrl() == "$url/acme/proj/hello"
        // the clone url ignores the configured server
        repo.getCloneUrl() == 'https://dev.azure.com/acme/proj/_git/hello'
        and:
        wireMock.allServeEvents.isEmpty()
    }

    @PendingFeature(reason = 'the clone url ignores the configured server')
    def 'should build the clone url from the configured server' () {
        expect:
        provider().getCloneUrl() == "$url/acme/proj/_git/hello"
    }

    def 'should read file content' () {
        given:
        final path = "$REPO_PATH/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&\$format=json&path=main.nf"
        stubJson(path, '''
            {"objectId": "f1", "gitObjectType": "blob", "commitId": "c1", "path": "/main.nf",
             "content": "println \\"Hello world!\\"\\n"}
            ''')

        expect:
        provider().readText('main.nf') == 'println "Hello world!"\n'
    }

    def 'should read file content at commit revision' () {
        given:
        final path = "$REPO_PATH/items?download=false&includeContent=true&includeContentMetadata=false&api-version=6.0&\$format=json&path=main.nf" +
            "&versionDescriptor.version=$COMMIT&versionDescriptor.versionType=commit"
        stubJson(path, '{"path": "/main.nf", "content": "at commit"}')

        expect:
        provider().setRevision(COMMIT).readText('main.nf') == 'at commit'
    }

    def 'should read binary content at branch revision' () {
        given:
        final path = "$REPO_PATH/items?download=true&includeContent=false&includeContentMetadata=false&api-version=6.0&path=img/logo.png" +
            "&versionDescriptor.version=feature%2Fx"
        stub(path, ok().withBody(BINARY))

        expect:
        provider().setRevision('feature/x').readBytes('img/logo.png') == BINARY
    }

    def 'should list branches following the continuation token' () {
        given:
        final path = "$REPO_PATH/refs?filter=heads&api-version=6.0"
        stub(path, okJson('''
            {"count": 1, "value": [{"name": "refs/heads/main", "objectId": "aaa"}]}
            ''').withHeader('x-ms-continuationtoken', 'next-1'))
        stubJson("$path&continuationToken=next-1", '''
            {"count": 1, "value": [{"name": "refs/heads/feature/x", "objectId": "bbb"}]}
            ''')

        when:
        def result = provider().getBranches()

        then:
        result == [
            new RepositoryProvider.BranchInfo('main', 'aaa'),
            new RepositoryProvider.BranchInfo('feature/x', 'bbb') ]
        wireMock.verify(2, getRequestedFor(urlPathEqualTo("$REPO_PATH/refs")))
    }

    def 'should list tags with the tag object of annotated tags' () {
        given:
        stubTags()

        expect:
        provider().getTags() == [
            new RepositoryProvider.TagInfo('v1.0', 't10'),
            new RepositoryProvider.TagInfo('v1.1', TAG_OBJECT) ]
    }

    @PendingFeature(reason = 'annotated tags resolve to the tag object instead of the commit')
    def 'should list tags with the commit of annotated tags' () {
        given:
        stubTags()

        expect:
        provider().getTags() == [
            new RepositoryProvider.TagInfo('v1.0', 't10'),
            new RepositoryProvider.TagInfo('v1.1', TAG_COMMIT) ]
    }

    def 'should list root directory' () {
        given:
        // the items list does not include the size
        stubJson("$REPO_PATH/items?recursionLevel=OneLevel&api-version=6.0&\$format=json", '''
            {"count": 3, "value": [
                {"objectId": "r0", "gitObjectType": "tree", "path": "/", "isFolder": true},
                {"objectId": "f1", "gitObjectType": "blob", "path": "/main.nf"},
                {"objectId": "d1", "gitObjectType": "tree", "path": "/docs", "isFolder": true}
            ]}
            ''')

        when:
        def result = provider().listDirectory('/', 1)

        then:
        result == [
            new RepositoryProvider.RepositoryEntry('docs', '/docs', RepositoryProvider.EntryType.DIRECTORY, 'd1', null),
            new RepositoryProvider.RepositoryEntry('main.nf', '/main.nf', RepositoryProvider.EntryType.FILE, 'f1', null) ]
    }

    def 'should list subdirectory recursively at revision' () {
        given:
        stubJson("$REPO_PATH/items?recursionLevel=Full&api-version=6.0&\$format=json&scopePath=docs&versionDescriptor.version=dev", '''
            {"count": 4, "value": [
                {"objectId": "d1", "gitObjectType": "tree", "path": "/docs", "isFolder": true},
                {"objectId": "f2", "gitObjectType": "blob", "path": "/docs/index.md"},
                {"objectId": "d2", "gitObjectType": "tree", "path": "/docs/guide", "isFolder": true},
                {"objectId": "f3", "gitObjectType": "blob", "path": "/docs/guide/intro.md"}
            ]}
            ''')

        when:
        def result = provider().setRevision('dev').listDirectory('docs', 2)

        then:
        result*.path == ['/docs/guide', '/docs/index.md', '/docs/guide/intro.md']
    }

    def 'should return an empty listing on errors' () {
        given:
        wireMock.stubFor(withAuth(get(urlPathEqualTo("$REPO_PATH/items")))
            .willReturn(aResponse().withStatus(404).withBody('{"message": "not found"}')))

        expect:
        provider().listDirectory('missing', 1) == []
        wireMock.verify(1, getRequestedFor(urlPathEqualTo("$REPO_PATH/items")))
    }

    def 'should report missing file' () {
        given:
        wireMock.stubFor(withAuth(get(urlPathEqualTo("$REPO_PATH/items")))
            .willReturn(aResponse().withStatus(404).withBody('{"message": "TF401174: The item could not be found"}')))

        when:
        provider().readText('missing.nf')

        then:
        thrown(AbortOperationException)
    }
}
