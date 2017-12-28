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
 * Created by mchatzou on 8/14/14.
 */
class BitBucketRepositoryProviderTest extends Specification {

    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def testBitbucketCloneURL() {

        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)

        when:
        def url = new BitbucketRepositoryProvider('pditommaso/tutorial',config).getCloneUrl()
        then:
        url == "https://${config.user}@bitbucket.org/pditommaso/tutorial.git".toString()
    }


    def testGetHomePage() {
        expect:
        new BitbucketRepositoryProvider('pditommaso/tutorial').getRepositoryUrl() == "https://bitbucket.org/pditommaso/tutorial"
    }


    @Requires( { System.getenv('NXF_BITBUCKET_ACCESS_TOKEN') } )
    def testReadContent() {

        given:
        def token = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')
        def config = new ProviderConfig('bitbucket').setAuth(token)

        when:
        def repo = new BitbucketRepositoryProvider('pditommaso/tutorial', config)
        def result = repo.readText('main.nf')

        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }
}
