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

package nextflow.processor

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
class ArchitectureTest extends Specification {

    @Unroll
    def 'should define the CPU architecture' () {

        when:
        def arch = new Architecture(VALUE)
        then:
        arch.dockerArch == DOCK
        arch.spackArch == SPACK

        where:
        VALUE                | DOCK             | SPACK
        'x86_64'             | 'linux/amd64'    | 'x86_64'
        'linux/x86_64'       | 'linux/amd64'    | 'x86_64'
        'amd64'              | 'linux/amd64'    | 'x86_64'
        'aarch64'            | 'linux/arm64'    | 'aarch64'
        'arm64'              | 'linux/arm64'    | 'aarch64'
        'linux/arm64/v8'     | 'linux/arm64'    | 'aarch64'
        'linux/arm64/v7'     | 'linux/arm64/v7' | null
    }

    def 'should define arch with map' () {
        when:
        def arch = new Architecture(VALUE)
        then:
        arch.dockerArch == DOCK
        arch.spackArch == SPACK

        where:
        VALUE                                  | DOCK             | SPACK
        [name: 'amd64', target: 'zen3']        | 'linux/amd64'    | 'zen3'
        [name: 'arm64', target: 'zen3']        | 'linux/arm64'    | 'zen3'
        [name: 'linux/x86_64', target: 'zen3'] | 'linux/amd64'    | 'zen3'
    }

    @Unroll
    def 'should normalize arch from name' () {
        expect:
        Architecture.ArchEntry.normalize(NAME) == EXPECTED

        where:
        NAME             | EXPECTED
        'x86_64'         | 'x86_64'
        'linux/x86_64'   | 'x86_64'
        'linux/arm64/v8' | 'arm64/v8'
        'linux/arm64/v7' | 'arm64/v7'
        'arm64'          | 'arm64'
        'aarch64'        | 'aarch64'
    }

    def 'should return dockerArch as string representation' () {
        expect:
        new Architecture('linux/x86_64').toString() == 'linux/amd64'
        new Architecture('arm64').toString() == 'linux/arm64'
        new Architecture('linux/amd64,linux/arm64').toString() == 'linux/amd64,linux/arm64'
    }

    def 'should parse multi-arch value' () {
        when:
        def arch = new Architecture('linux/amd64,linux/arm64')
        then:
        arch.dockerArch == 'linux/amd64'
        arch.spackArch == 'x86_64'
        arch.platforms() == ['linux/amd64', 'linux/arm64']
    }

    def 'should parse multi-arch with spaces' () {
        when:
        def arch = new Architecture('linux/amd64, linux/arm64')
        then:
        arch.platforms() == ['linux/amd64', 'linux/arm64']
    }

    def 'should return containerPlatform as comma-separated string' () {
        expect:
        new Architecture('linux/amd64').containerPlatform() == 'linux/amd64'
        new Architecture('linux/amd64,linux/arm64').containerPlatform() == 'linux/amd64,linux/arm64'
    }

    def 'should define multi-arch with map' () {
        when:
        def arch = new Architecture([name: 'amd64,arm64', target: 'zen3'])
        then:
        arch.dockerArch == 'linux/amd64'
        arch.spackArch == 'zen3'
        arch.platforms() == ['linux/amd64', 'linux/arm64']
        arch.containerPlatform() == 'linux/amd64,linux/arm64'
    }

    def 'should use primary arch for dockerArch and spackArch with multi-arch' () {
        when:
        def arch = new Architecture('arm64,amd64')
        then:
        arch.dockerArch == 'linux/arm64'
        arch.spackArch == 'aarch64'
        arch.platforms() == ['linux/arm64', 'linux/amd64']
    }
}
