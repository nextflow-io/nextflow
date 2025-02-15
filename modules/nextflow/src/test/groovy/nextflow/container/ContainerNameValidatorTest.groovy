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

package nextflow.container

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerNameValidatorTest extends Specification {

    def 'should match image name' () {

        expect:
        ContainerNameValidator.isValidImageName(IMG) == EXPECTED

        where:
        EXPECTED    | IMG
        true        | 'image-name'
        true        | 'image-name:latest'
        true        | 'image-name:123'
        true        | 'registry.example.com/org/image-name'
        true        | 'registry/org/image-name'
        true        | 'registry/image-name'
        true        | 'image-name'
        true        | 'registry.example.com/org/image-name:version'
        true        | 'registry.example.com/org/image-name@sha256:c0e9560cda118f9ec63ddefb4a173a2b2a0347082d7dff7dc14272e7841a5b5a'
        true        | 'registry/org/image-name:version'
        true        | 'registry/org/image-name@sha256:c0e9560cda118f9ec63ddefb4a173a2b2a0347082d7dff7dc14272e7841a5b5a'
        true        | 'registry/image-name:version'
        true        | 'registry/image-name@sha256:c0e9560cda118f9ec63ddefb4a173a2b2a0347082d7dff7dc14272e7841a5b5a'
        true        | 'registry:8000/image-name:version'
        true        | 'registry:8000/image-name@sha256:c0e9560cda118f9ec63ddefb4a173a2b2a0347082d7dff7dc14272e7841a5b5a'
        true        | 'image-name:version'
        true        | 'image-name:VERSION'
        true        | 'image-name@sha256:c0e9560cda118f9ec63ddefb4a173a2b2a0347082d7dff7dc14272e7841a5b5a'
        true        | 'wave.seqera.io/wt/b3f144027b00/biocontainers/control-freec:11.6b--hdbdd923_0'
        and:
        false       | 'image-name:version '
        false       | 'IMAGE-NAME:123'
        false       | '/some/image'
        false       | 'http://registry/image-name:version'
        false       | 'registry:/image-name:version'

    }
}
