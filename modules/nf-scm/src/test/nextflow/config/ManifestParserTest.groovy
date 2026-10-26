/*
 * Copyright 2013-2026, Seqera Labs
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

import spock.lang.Specification

class ManifestParserTest extends Specification {

    def 'should read literal manifest values from a block' () {
        when:
        def manifest = ManifestParser.parse('''
            manifest {
                description = 'A pipeline'
                mainScript = 'other.nf'
                defaultBranch = 'dev'
                recurseSubmodules = true
                gitmodules = ['a', 'b']
                contributors = [ [name: 'Alice', contribution: ['author']] ]
            }
            ''')
        then:
        manifest.description == 'A pipeline'
        manifest.mainScript == 'other.nf'
        manifest.defaultBranch == 'dev'
        manifest.recurseSubmodules
        manifest.gitmodules == ['a', 'b']
        manifest.contributors*.name == ['Alice']
    }

    def 'should read dotted manifest assignments' () {
        when:
        def manifest = ManifestParser.parse('''
            manifest.homePage = 'https://example.org'
            process.cpus = 2
            ''')
        then:
        manifest.homePage == 'https://example.org'
    }

    def 'should ignore non-literal values and other scopes' () {
        when:
        def manifest = ManifestParser.parse('''
            params.suffix = 'dev'
            manifest {
                version = "1.0-${params.suffix}"
                mainScript = 'main2.nf'
            }
            profiles {
                test { manifest.mainScript = 'test.nf' }
            }
            ''')
        then:
        manifest.version == null
        manifest.mainScript == 'main2.nf'
    }

    def 'should not evaluate the config' () {
        when:
        def manifest = ManifestParser.parse('''
            manifest.description = 'safe'
            params.x = { throw new IllegalStateException('evaluated') }()
            ''')
        then:
        manifest.description == 'safe'
    }

    def 'should fail on invalid syntax' () {
        when:
        ManifestParser.parse('manifest { description = ')
        then:
        thrown(IllegalArgumentException)
    }
}
