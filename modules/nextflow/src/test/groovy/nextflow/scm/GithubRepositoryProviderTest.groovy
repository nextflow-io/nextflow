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

class GithubRepositoryProviderTest extends Specification {

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testGitCloneUrl() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)

        when:
        def url = new GithubRepositoryProvider('nextflow-io/hello',config).getCloneUrl()
        then:
        url == 'https://github.com/nextflow-io/hello.git'

    }

    def testGetHomePage() {
        expect:
        new GithubRepositoryProvider('nextflow-io/hello').getRepositoryUrl() == "https://github.com/nextflow-io/hello"
    }

    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testReadContent() {

        given:
        def token = System.getenv('NXF_GITHUB_ACCESS_TOKEN')
        def config = new ProviderConfig('github').setAuth(token)

        when:
        def repo = new GithubRepositoryProvider('nextflow-io/hello', config)
        def result = repo.readText('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }
}

