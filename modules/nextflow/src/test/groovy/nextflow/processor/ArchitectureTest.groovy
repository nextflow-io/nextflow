/*
 * Copyright 2023, Pawsey Supercomputing Research Centre
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
        arch.platform == PLAT
        arch.arch == ARCH
        arch.target == TAR
        arch.dockerArch == DOCK
        arch.spackArch == SPACK

        where:
        VALUE                                  | PLAT     | ARCH        | TAR    | DOCK             | SPACK
        'x86_64'                               | null     | 'x86_64'    | null   | 'linux/amd64'    | 'x86_64'
        'linux/x86_64'                         | 'linux'  | 'x86_64'    | null   | 'linux/amd64'    | 'x86_64'
        'amd64'                                | null     | 'amd64'     | null   | 'linux/amd64'    | 'x86_64'
        'aarch64'                              | null     | 'aarch64'   | null   | 'linux/arm64'    | 'aarch64'
        'arm64'                                | null     | 'arm64'     | null   | 'linux/arm64'    | 'aarch64'
        'linux/arm64/v8'                       | 'linux'  | 'arm64/v8'  | null   | 'linux/arm64'    | 'aarch64'
        'linux/arm64/v7'                       | 'linux'  | 'arm64/v7'  | null   | 'linux/arm64/v7' | null
    }

    def 'should define arch with map' () {
        when:
        def arch = new Architecture(VALUE)
        then:
        arch.platform == PLAT
        arch.arch == ARCH
        arch.target == TAR
        arch.dockerArch == DOCK
        arch.spackArch == SPACK

        where:
        VALUE                                  | PLAT     | ARCH        | TAR    | DOCK             | SPACK
        [name: 'amd64', target: 'zen3']        | null     | 'amd64'     | 'zen3' | 'linux/amd64'    | 'zen3'
        [name: 'arm64', target: 'zen3']        | null     | 'arm64'     | 'zen3' | 'linux/arm64'    | 'zen3'
        [name: 'linux/x86_64', target: 'zen3'] | 'linux'  | 'x86_64'    | 'zen3' | 'linux/amd64'    | 'zen3'
    }
}
