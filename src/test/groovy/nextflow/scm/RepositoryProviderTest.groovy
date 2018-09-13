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

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class RepositoryProviderTest extends Specification {

    def 'should create repository provider object' () {

        def provider

        when:
        provider = RepositoryProvider.create(new ProviderConfig('github'),'project/x')
        then:
        provider instanceof GithubRepositoryProvider
        provider.endpointUrl == 'https://api.github.com/repos/project/x'

        when:
        provider = RepositoryProvider.create(new ProviderConfig('gitlab'),'project/y')
        then:
        provider instanceof GitlabRepositoryProvider
        provider.endpointUrl == 'https://gitlab.com/api/v4/projects/project%2Fy'

        when:
        provider = RepositoryProvider.create(new ProviderConfig('bitbucket'),'project/z')
        then:
        provider instanceof BitbucketRepositoryProvider
        provider.endpointUrl == 'https://bitbucket.org/api/2.0/repositories/project/z'

        when:
        provider = RepositoryProvider.create(new ProviderConfig('local', [path:'/user/data']),'local/w')
        then:
        provider.endpointUrl == 'file:/user/data/w'
    }

    def 'should set credentials' () {

        given:
        def config = Mock(ProviderConfig)
        def provider = Spy(RepositoryProvider)
        provider.config = config

        when:
        provider.setCredentials('pditommaso', 'secret1')
        then:
        1 * config.setUser('pditommaso')
        1 * config.setPassword('secret1')

    }
}
