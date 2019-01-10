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

package nextflow.config

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ManifestTest extends Specification {

    def 'should check manifest object' () {

        given:
        def MAN = [author: 'pablo', nextflowVersion: '1.2.3', name: 'foo']

        when:
        def manifest = new Manifest(MAN)
        then:
        manifest.with {
            author == 'pablo'
            nextflowVersion == '1.2.3'
            name == 'foo'
        }

    }

    def 'should check empty manifest' () {

        // check empty manifest
        when:
        def manifest = new Manifest(new ConfigObject())
        then:
        manifest.with {
            homePage == null
            defaultBranch == 'master'
            description == null
            author == null
            mainScript == 'main.nf'
            gitmodules == null
            nextflowVersion == null
            version == null
            name == null
        }

    }


}
