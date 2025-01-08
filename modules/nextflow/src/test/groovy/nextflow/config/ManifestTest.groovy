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

package nextflow.config

import nextflow.exception.AbortOperationException
import spock.lang.Specification

import static nextflow.config.Manifest.ContributionType
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ManifestTest extends Specification {

    def 'should check manifest object' () {

        given:
        def MAP = [
            author: 'pablo',
            contributors: [
                [
                    name: 'Alice', 
                    affiliation: 'University', 
                    email: 'alice@university.edu',
                    contribution: ['author', 'maintainer'],
                    orcid: 'https://orcid.org/0000-0000-0000-0000'
                ],
                [
                    name: 'Bob', 
                    affiliation: 'Company', 
                    email: 'bob@company.com',
                    contribution: ['contributor'],
                ]
            ],
            nextflowVersion: '1.2.3',
            name: 'foo',
            organization: 'My Organization',
            icon: 'icon.png',
            docsUrl: 'https://docs.io',
            license: 'Apache v2'
        ]
        when:
        def manifest = new Manifest(MAP)
        then:
        manifest.author == 'pablo'
        manifest.contributors.size() == 2
        manifest.contributors[0].name == 'Alice'
        manifest.contributors[0].affiliation == 'University'
        manifest.contributors[0].email == 'alice@university.edu'
        manifest.contributors[0].contribution == [ContributionType.AUTHOR, ContributionType.MAINTAINER] as Set
        manifest.contributors[0].orcid == 'https://orcid.org/0000-0000-0000-0000'
        manifest.contributors[1].name == 'Bob'
        manifest.contributors[1].affiliation == 'Company'
        manifest.contributors[1].email == 'bob@company.com'
        manifest.contributors[1].contribution == [ContributionType.CONTRIBUTOR] as Set
        manifest.nextflowVersion == '1.2.3'
        manifest.name == 'foo'
        manifest.organization == 'My Organization'
        manifest.icon == 'icon.png'
        manifest.docsUrl == 'https://docs.io'
        manifest.license == 'Apache v2'

    }

    def 'should check empty manifest' () {

        when:
        def manifest = new Manifest(new ConfigObject())
        then:
        manifest.homePage == null
        manifest.defaultBranch == null
        manifest.description == null
        manifest.author == null
        manifest.contributors == []
        manifest.mainScript == 'main.nf'
        manifest.gitmodules == null
        manifest.nextflowVersion == null
        manifest.version == null
        manifest.name == null
        manifest.docsUrl == null
        manifest.organization == null
        manifest.icon == null
        manifest.license == null

    }

    def 'should convert manifest to map' () {

        when:
        def MAP = [
            name: 'Alice', 
            affiliation: 'University', 
            email: 'alice@university.edu',
            contribution: ['author', 'maintainer'],
            orcid: 'https://orcid.org/0000-0000-0000-0000'
        ]
        then:
        new Manifest.Contributor(MAP).toMap() == [
            name: 'Alice', 
            affiliation: 'University', 
            email: 'alice@university.edu',
            github: null,
            contribution: ['author', 'maintainer'],
            orcid: 'https://orcid.org/0000-0000-0000-0000'
        ]
    }

    def 'should throw error on invalid manifest' () {
        when:
        def manifest = new Manifest([
            contributors: [ 'Alice' ]
        ])
        manifest.contributors
        then:
        thrown(AbortOperationException)

        when:
        manifest = new Manifest([
            contributors: [[
                name: 'Alice',
                contribution: 'author'
            ]]
        ])
        manifest.contributors
        then:
        thrown(AbortOperationException)

        when:
        manifest = new Manifest([
            contributors: [[
                name: 'Alice',
                contribution: [ 'owner' ]
            ]]
        ])
        manifest.contributors
        then:
        thrown(AbortOperationException)
    }

}
