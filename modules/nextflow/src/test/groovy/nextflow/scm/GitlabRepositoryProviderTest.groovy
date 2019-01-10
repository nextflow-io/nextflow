/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

}
