/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.scm

import spock.lang.Requires
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
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
