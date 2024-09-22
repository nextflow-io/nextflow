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
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ManifestTest extends Specification {

    def 'should check manifest object' () {

        given:
        def MAN = [
            author: 'pablo',
            contributors: [
                [
                    name: 'Alice', 
                    affiliation: 'University', 
                    email: 'alice@university.edu',
                    contribution: ['author', 'maintainer'],
                ],
                [
                    name: 'Bob', 
                    affiliation: 'Company', 
                    email: 'bob@company.com',
                    contribution: ['maintainer'],
                ]
            ],
            nextflowVersion: '1.2.3',
            name: 'foo',
            organisation: 'My Organisation',
            icon: 'icon.png',
            docsUrl: 'https://docs.io',
            license: 'Apache v2'
        ]
        when:
        def manifest = new Manifest(MAN)
        then:
        manifest.author == 'pablo'
        manifest.contributors == [
            new Manifest.Contributor([
                name: 'Alice', 
                affiliation: 'University', 
                email: 'alice@university.edu',
                contribution: ['author', 'maintainer'],
            ]),
            new Manifest.Contributor([
                name: 'Bob', 
                affiliation: 'Company', 
                email: 'bob@company.com',
                contribution: ['maintainer'],
            ])
        ]
        manifest.nextflowVersion == '1.2.3'
        manifest.name == 'foo'
        manifest.organisation == 'My Organisation'
        manifest.icon == 'icon.png'
        manifest.docsUrl == 'https://docs.io'
        manifest.license == 'Apache v2'

    }

    def 'should check empty manifest' () {

        when:
        def manifest = new Manifest(new ConfigObject())
        then:
        manifest.homePage == null
        manifest.defaultBranch == 'master'
        manifest.description == null
        manifest.author == null
        manifest.contributors == []
        manifest.mainScript == 'main.nf'
        manifest.gitmodules == null
        manifest.nextflowVersion == null
        manifest.version == null
        manifest.name == null
        manifest.docsUrl == null
        manifest.organisation == null
        manifest.icon == null
        manifest.license == null

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
