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
        true        | 'registry/org/image-name:version'
        true        | 'registry/image-name:version'
        true        | 'registry:8000/image-name:version'
        true        | 'image-name:version'
        true        | 'image-name:VERSION'
        and:
        false       | 'IMAGE-NAME:123'
        false       | '/some/image'
        false       | 'http://registry/image-name:version'
        false       | 'registry:/image-name:version'

    }
}
