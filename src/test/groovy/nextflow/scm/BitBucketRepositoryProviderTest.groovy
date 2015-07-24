/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
        def (user,pwd) = System.getenv('NXF_BITBUCKET_ACCESS_TOKEN')?.tokenize(':')

        when:
        def url = new BitbucketRepositoryProvider(pipeline: 'pditommaso/tutorial', user: user, pwd: pwd).getCloneUrl()
        then:
        url == "https://${user}@bitbucket.org/pditommaso/tutorial.git".toString()
    }


    def testGetHomePage() {
        expect:
        new BitbucketRepositoryProvider(pipeline: 'pditommaso/tutorial').getHomePage() == "https://bitbucket.org/pditommaso/tutorial"
    }


    def testReadContent() {

        when:
        def repo = new BitbucketRepositoryProvider(pipeline: 'pditommaso/tutorial')
        def result = repo.readText('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }
}
