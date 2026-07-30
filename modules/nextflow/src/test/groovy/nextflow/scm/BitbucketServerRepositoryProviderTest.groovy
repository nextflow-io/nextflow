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

import org.eclipse.jgit.transport.CredentialItem
import org.eclipse.jgit.transport.URIish
import spock.lang.Ignore
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Piotr Faba <piotr.faba@ardigen.com>
 */
class BitbucketServerRepositoryProviderTest extends Specification {

    static final String CONFIG = '''
        providers {

            bbserver {
                server = 'https://bitbucket.server.com'
                endpoint = 'https://bitbucket.server.com'
                platform = 'bitbucketserver'
                user = 'myname'
                password = 'mypassword'
            }

        }
        '''

    @Unroll
    def 'should return endpoint url' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider(NAME, obj).getEndpointUrl() == ENDPOINT

        where:
        NAME                                    | ENDPOINT
        'pditommaso/hello'                      | 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello'
        'scm/pditommaso/hello'                  | 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello'
        'projects/DA/repos/crispr-pipeline'     | 'https://bitbucket.server.com/rest/api/1.0/projects/DA/repos/crispr-pipeline'
    }

    @Unroll
    def 'should return project URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider(NAME, obj).getRepositoryUrl() == REPO_URL

        where:
        NAME                                | REPO_URL
//        'pditommaso/hello'                  | 'https://bitbucket.server.com/pditommaso/hello'
        '/pditommaso/hello'                 | 'https://bitbucket.server.com/pditommaso/hello'
//        and:
//        'scm/pditommaso/hello'              | 'https://bitbucket.server.com/scm/pditommaso/hello'
//        '/scm/pditommaso/hello'             | 'https://bitbucket.server.com/scm/pditommaso/hello'
//        and:
//        'projects/DA/repos/crispr-pipeline' | 'https://bitbucket.server.com/projects/DA/repos/crispr-pipeline'
//        '/projects/DA/repos/crispr-pipeline'| 'https://bitbucket.server.com/projects/DA/repos/crispr-pipeline'
    }

    def 'should return content URL' () {

        given:
        def config = new ConfigSlurper().parse(CONFIG)
        def obj = new ProviderConfig('bitbucketserver', config.providers.bbserver as ConfigObject)

        expect:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf'

        and:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj)
                .setRevision('foo')
                .getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf?at=foo'
        and:
        new BitbucketServerRepositoryProvider('pditommaso/hello', obj)
                .setRevision('test/branch+with&strangecharacters')
                .getContentUrl('main.nf') == 'https://bitbucket.server.com/rest/api/1.0/projects/pditommaso/repos/hello/raw/main.nf?at=test%2Fbranch%2Bwith%26strangecharacters'

    }

    @Ignore
    @Requires( { System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN') } )
    def 'should list branches' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN')
        def items = token.tokenize(':')

        def config = items.size() > 1
            ? new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketserver']).setUser(items[0]).setToken(items[1])
            : new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketserver']).setToken(token)
        and:
        def repo = new BitbucketServerRepositoryProvider('scm/hello/hello', config)

        when:
        def result = repo.getBranches()
        then:
        result.contains( new RepositoryProvider.BranchInfo('master', 'c62df3d9c2464adcaa0fb6c978c8e32e2672b191') )
    }

    @Ignore
    @Requires( { System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN') } )
    def 'should list tags' () {
        given:
        def token = System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN')
        def config = new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketserver']).setAuth(token)
        and:
        def repo = new BitbucketServerRepositoryProvider('scm/hello/hello', config)

        when:
        def result = repo.getTags()
        then:
        result.contains( new RepositoryProvider.TagInfo('v1.0', 'c62df3d9c2464adcaa0fb6c978c8e32e2672b191') )
    }

    @Unroll
    def 'should validate hasCredentials' () {
        given:
        def provider = new BitbucketServerRepositoryProvider('proj/repo', PROVIDER)

        expect:
        provider.hasCredentials() == EXPECTED

        where:
        EXPECTED | PROVIDER
        false    | new ProviderConfig('bbs', [platform:'bitbucketserver'])
        false    | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo')
        false    | new ProviderConfig('bbs', [platform:'bitbucketserver']).setPassword('bar')
        true     | new ProviderConfig('bbs', [platform:'bitbucketserver']).setToken('xyz')
        true     | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setPassword('bar')
        true     | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setToken('xyz')
        true     | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setPassword('bar').setToken('xyz')
    }

    @Unroll
    def 'should validate getAuth' () {
        given:
        def provider = new BitbucketServerRepositoryProvider('proj/repo', PROVIDER)

        expect:
        provider.getAuth() == EXPECTED as String[]

        where:
        EXPECTED                                                       | PROVIDER
        null                                                           | new ProviderConfig('bbs', [platform:'bitbucketserver'])
        ["Authorization", "Bearer xyz"]                                | new ProviderConfig('bbs', [platform:'bitbucketserver']).setToken('xyz')
        ["Authorization", "Basic ${"foo:bar".bytes.encodeBase64()}"]   | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setPassword('bar')
        ["Authorization", "Bearer xyz"]                                | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setToken('xyz')
        ["Authorization", "Bearer xyz"]                                | new ProviderConfig('bbs', [platform:'bitbucketserver']).setUser('foo').setPassword('bar').setToken('xyz')
    }

    def 'should pass token as password in getGitCredentials' () {
        given:
        def config = new ProviderConfig('bbs', [platform:'bitbucketserver'])
                .setUser('foo')
                .setToken('xyz')
        def provider = new BitbucketServerRepositoryProvider('proj/repo', config)
        def user = new CredentialItem.Username()
        def pass = new CredentialItem.Password()

        when:
        def creds = provider.getGitCredentials()
        creds.get(new URIish('https://bitbucket.server.com/scm/proj/repo.git'), user, pass)

        then:
        user.value == 'foo'
        new String(pass.value) == 'xyz'
    }

    def 'should use token-only in getGitCredentials when no user is set' () {
        given:
        def config = new ProviderConfig('bbs', [platform:'bitbucketserver'])
                .setAuth('xyz')
        def provider = new BitbucketServerRepositoryProvider('proj/repo', config)
        def user = new CredentialItem.Username()
        def pass = new CredentialItem.Password()

        when:
        def creds = provider.getGitCredentials()
        creds.get(new URIish('https://bitbucket.server.com/scm/proj/repo.git'), user, pass)

        then:
        user.value == ''
        new String(pass.value) == 'xyz'
    }

    def 'should fall back to password in getGitCredentials when token absent' () {
        given:
        def config = new ProviderConfig('bbs', [platform:'bitbucketserver'])
                .setUser('foo')
                .setPassword('bar')
        def provider = new BitbucketServerRepositoryProvider('proj/repo', config)
        def user = new CredentialItem.Username()
        def pass = new CredentialItem.Password()

        when:
        def creds = provider.getGitCredentials()
        creds.get(new URIish('https://bitbucket.server.com/scm/proj/repo.git'), user, pass)

        then:
        user.value == 'foo'
        new String(pass.value) == 'bar'
    }

    @Requires( { System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN') } )
    def 'should list root directory contents'() {
        given:
        def token = System.getenv('NXF_BITBUCKET_SERVER_ACCESS_TOKEN')
        def config = new ProviderConfig('bbs', [server:'http://slurm.seqera.io:7990', platform:'bitbucketserver']).setAuth(token)
        def repo = new BitbucketServerRepositoryProvider('scm/hello/hello', config)

        when:
        repo.listDirectory("/", 1)

        then:
        thrown UnsupportedOperationException
    }

}
