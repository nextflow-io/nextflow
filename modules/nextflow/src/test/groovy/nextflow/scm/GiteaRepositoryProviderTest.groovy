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

import spock.lang.IgnoreIf
import spock.lang.Requires
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Akira Sekiguchi <pachiras.yokohama@gmail.com>
 */
class GiteaRepositoryProviderTest extends Specification {

    def 'should return repo url' () {
        given:
        def obj = new ProviderConfig('gitea')

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getEndpointUrl() == 'https://gitea.com/api/v1/repos/pditommaso/hello'
    }

    def 'should return project URL' () {
        given:
        def obj = new ProviderConfig('gitea')

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj).getRepositoryUrl() == 'https://gitea.com/pditommaso/hello'
    }

    def 'should return content URL' () {
        given:
        def obj = new ProviderConfig('gitea')

        expect:
        new GiteaRepositoryProvider('pditommaso/hello', obj)
                .getContentUrl('main.nf') == 'https://gitea.com/api/v1/repos/pditommaso/hello/raw/main.nf'
        and:
        new GiteaRepositoryProvider('pditommaso/hello', obj)
                .setRevision('12345')
                .getContentUrl('main.nf') == 'https://gitea.com/api/v1/repos/pditommaso/hello/raw/main.nf?ref=12345'

    }

    @Unroll
    def 'should validate hasCredentials' () {
        given:
        def provider = new GiteaRepositoryProvider('pditommaso/tutorial', PROVIDER_CONFIG)

        expect:
        provider.hasCredentials() == EXPECTED

        where:
        EXPECTED    | PROVIDER_CONFIG
        false       | new ProviderConfig('gitea')
        false       | new ProviderConfig('gitea').setUser('foo')
        true        | new ProviderConfig('gitea').setUser('foo').setPassword('bar')
        true        | new ProviderConfig('gitea').setToken('xyz')
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_GITEA_ACCESS_TOKEN')})
    def 'should read file content'() {
        given:
        def token =  System.getenv('NXF_GITEA_ACCESS_TOKEN')
        def config = new ProviderConfig('gitea') .setAuth(token)

        when:
        def repo = new GiteaRepositoryProvider('pditommaso/test-hello', config)
        def result = repo.readText('README.md')
        then:
        result.contains('Basic Nextflow script')

//        when:
//        repo = new GiteaRepositoryProvider('test-org/nextflow-ci-repo', config)
//                        .setRevision('foo')
//        result = repo.readText('README.md')
//        then:
//        result.contains("foo branch")
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_GITEA_ACCESS_TOKEN')})
    def 'should read bytes gitea content'() {
        given:
        def token =  System.getenv('NXF_GITEA_ACCESS_TOKEN')
        def config = new ProviderConfig('gitea') .setAuth(token)
        def repo = new GiteaRepositoryProvider('pditommaso/test-hello', config)
        and:
        def DATA = this.class.getResourceAsStream('/test-asset.bin').bytes
        
        when:
        def result = repo.readBytes('test/test-asset.bin')

        then:
        result == DATA
    }

    @IgnoreIf({System.getenv('NXF_SMOKE')})
    @Requires({System.getenv('NXF_GITEA_ACCESS_TOKEN')})
    def 'should read bytes file content'() {
        given:
        def token =  System.getenv('NXF_GITEA_ACCESS_TOKEN')
        def config = new ProviderConfig('gitea').setAuth(token)

        when:
        def repo = new GiteaRepositoryProvider('pditommaso/test-hello', config)
        def result = repo.readBytes('docs/images/nf-core-rnaseq_logo_light.png')

        then:
        result.length == 22915
        result.sha256() == '7a396344498750f614155f6e4f38b7d6ca98ced45daf0921b64acf73b18efaf4'
    }
}
