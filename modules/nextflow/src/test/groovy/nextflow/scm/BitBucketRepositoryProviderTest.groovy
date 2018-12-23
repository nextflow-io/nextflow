/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.Requires
import spock.lang.Specification

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
