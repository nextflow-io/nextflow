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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerConfigTest extends Specification {

    @Unroll
    def 'should return env whitelist for=#VAL' () {
        when:
        def cfg = new ContainerConfig(envWhitelist: VAL)
        then:
        cfg.getEnvWhitelist() == EXPECTED

        where:
        VAL         | EXPECTED
        null        | []
        ''          | []
        'FOO'       | ['FOO']
        'FOO,BAR'   | ['FOO','BAR']
        'A ,, B,C ' | ['A','B','C']
        ['X','Y']   | ['X','Y']

    }

    def 'should validate legacy entry point' () {

        when:
        def cfg = new ContainerConfig(OPTS, ENV)
        then:
        cfg.entrypointOverride() == EXPECTED
        
        where:
        OPTS                            | ENV          | EXPECTED
        [:]                             | [:]          | false
        [entrypointOverride: false]     | [:]          | false
        [entrypointOverride: true]      | [:]          | true
        and:
        [:]                             | [NXF_CONTAINER_ENTRYPOINT_OVERRIDE: 'true']  | true
        [entrypointOverride: false]     | [NXF_CONTAINER_ENTRYPOINT_OVERRIDE: 'true']  | false

    }

    @Unroll
    def 'should validate oci mode and direct mode' () {

        when:
        def cfg = new ContainerConfig(OPTS)
        then:
        cfg.isSingularityOciMode() == OCI_MODE
        cfg.canRunOciImage() == AUTO_PULL

        where:
        OPTS                                        | OCI_MODE  | AUTO_PULL
        [:]                                         | false     | false
        [oci:true]                                  | false     | false
        [oci:false]                                 | false     | false
        [ociMode:true]                              | false     | false
        and:
        [engine:'docker', oci:true]                 | false     | false
        [engine:'singularity']                      | false     | false
        [engine:'singularity', oci:false]           | false     | false
        [engine:'singularity', ociAutoPull:false]   | false     | false
        and:
        [engine:'singularity', oci:true]            | true      | true
        [engine:'singularity', ociMode:true]        | true      | true
        [engine:'apptainer', oci:true]              | false     | false
        [engine:'apptainer', ociMode:true]          | false     | false
        and:
        [engine:'singularity', ociAutoPull:true]    | false         | true
        [engine:'apptainer', ociAutoPull:true]      | false         | true

    }

    def 'should get fusion options' () {
        when:
        def cfg = new ContainerConfig(OPTS)

        then:
        cfg.fusionOptions() == EXPECTED
        
        where:
        OPTS                                            | EXPECTED
        [:]                                             | null
        [engine:'docker']                               | '--rm --privileged'
        [engine:'podman']                               | '--rm --privileged'
        and:
        [engine: 'singularity']                         | null
        [engine: 'singularity', ociMode:true]           | '-B /dev/fuse'
        [engine: 'singularity', ociAutoPull: true]      | null
        [engine: 'apptainer', oci:true]                 | null
        and:
        [engine:'docker', fusionOptions:'--cap-add foo']| '--cap-add foo'
        [engine:'podman', fusionOptions:'--cap-add bar']| '--cap-add bar'
        and:
        [engine:'sarus', fusionOptions:'--other']       | '--other'
    }

}
