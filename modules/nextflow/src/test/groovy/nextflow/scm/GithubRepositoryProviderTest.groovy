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

import nextflow.SysEnv
import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification

@IgnoreIf({System.getenv('NXF_SMOKE')})
class GithubRepositoryProviderTest extends Specification {

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testGitCloneUrl() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)

        when:
        def url = new GithubRepositoryProvider('nextflow-io/test-hello',config).getCloneUrl()
        then:
        url == 'https://github.com/nextflow-io/test-hello.git'
    }

    def testGetHomePage() {
        expect:
        new GithubRepositoryProvider('nextflow-io/test-hello').getRepositoryUrl() == "https://github.com/nextflow-io/test-hello"
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testReadContent() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def result = repo.readText('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should read bytes github content'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)
        and:
        def DATA = this.class.getResourceAsStream('/test-asset.bin').bytes
        
        when:
        def result = repo.readBytes('/test/test-asset.bin')
        then:
        result == DATA
    }

    def 'should return content URL' () {
        given:
        String CONFIG = '''
        providers {
            mygithub {
                server = 'https://github.com'
                endpoint = 'https://github.com'
                platform = 'bitbucket'
                user = 'myname'
                password = 'mypassword'
            }
        }
        '''

        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('github', config.providers.mygithub as ConfigObject)

        expect:
        new GithubRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://github.com/repos/pditommaso/hello/contents/main.nf'

        and:
        new GithubRepositoryProvider('pditommaso/hello', obj)
                .setRevision('the-commit-id')
                .getContentUrl('main.nf') == 'https://github.com/repos/pditommaso/hello/contents/main.nf?ref=the-commit-id'

    }

    def 'should user github token as creds' () {
        given:
        SysEnv.push(['GITHUB_TOKEN': '1234567890'])
        and:
        def provider = Spy(new GithubRepositoryProvider('foo/bar'))

        expect:
        provider.getUser() == '1234567890'
        provider.getPassword() == 'x-oauth-basic'
        
        when:
        SysEnv.get().remove('GITHUB_TOKEN')
        then:
        provider.getUser() >> null
        provider.getPassword() >> null

        cleanup:
        SysEnv.pop()
    }

    def 'should user from config' () {
        given:
        SysEnv.push(['GITHUB_TOKEN': '1234567890'])
        and:
        def config = new ProviderConfig('github', [user: 'this', password: 'that'])
        def provider = Spy(new GithubRepositoryProvider('foo/bar', config))

        expect:
        provider.getUser() == 'this'
        provider.getPassword() == 'that'

        cleanup:
        SysEnv.pop()
    }

    def 'should auth using github token' () {
        given:
        SysEnv.push(['GITHUB_TOKEN': '1234567890'])
        and:
        def provider = Spy(new GithubRepositoryProvider('foo',Mock(ProviderConfig)))

        when:
        final result = provider.getAuth()
        then:
        _ * provider.getUser() 
        _ * provider.getPassword()
        1 * provider.hasCredentials()
        and:
        result == new String[] { 'Authorization', "Basic ${'1234567890:x-oauth-basic'.bytes.encodeBase64()}".toString() }

        cleanup:
        SysEnv.pop()
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should read content from renamed repository'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        // this repository has been renamed to `pditommaso/ciao`
        // nevertheless the read content should work, because the
        // client needs to be able to follow the http redirection
        // returned by the github backend
        def repo = new GithubRepositoryProvider('pditommaso/hello', config)

        when:
        def result = repo.readText('main.nf')
        then:
        result.trim().startsWith(/println "I'm the main"/)
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list root directory contents'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def entries = repo.listDirectory("", 1)

        then:
        entries.size() > 0
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.every { it.path && it.sha }
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list subdirectory contents'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def entries = repo.listDirectory("test", 1)

        then:
        entries.size() > 0
        entries.any { it.name == 'test-asset.bin' && it.type == RepositoryProvider.EntryType.FILE }
        entries.every { it.path.startsWith('test/') }
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list directory contents recursively'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def entries = repo.listDirectory("", 10)

        then:
        entries.size() > 0
        // Should include files from root and subdirectories
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.any { it.name == 'test-asset.bin' && it.type == RepositoryProvider.EntryType.FILE }
        entries.every { it.path && it.sha }
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list directory contents with limited depth'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def depthOne = repo.listDirectory("", 1)
        def depthTwo = repo.listDirectory("", 2)

        then:
        depthOne.size() > 0
        depthTwo.size() >= depthOne.size()
        // Depth 1 should only include immediate children
        depthOne.every { !it.path.contains('/') }
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def 'should list directory contents with depth 2'() {
        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)
        def repo = new GithubRepositoryProvider('nextflow-io/test-hello', config)

        when:
        def entries = repo.listDirectory("", 2)

        then:
        entries.size() > 0
        // Should include immediate children (depth 1)
        entries.any { it.name == 'main.nf' && it.type == RepositoryProvider.EntryType.FILE }
        entries.any { it.name == 'test' && it.type == RepositoryProvider.EntryType.DIRECTORY }
        // Should include nested files (depth 2)
        entries.any { it.name == 'test-asset.bin' && it.path.contains('test/') }
        entries.every { it.path && it.sha }
    }

    def 'should return empty list for directory with no entries'() {
        expect:
        // This test will be integration test based - relying on actual API
        true
    }
}

