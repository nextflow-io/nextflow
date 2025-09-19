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

import nextflow.SysEnv
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
        def cfg = new DockerConfig(envWhitelist: VAL)
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
        SysEnv.push(ENV)
        def cfg = new DockerConfig(OPTS)
        def result = cfg.entrypointOverride()
        SysEnv.pop()
        then:
        result == EXPECTED

        where:
        OPTS    | ENV          | EXPECTED
        [:]     | [:]          | false
        and:
        [:]     | [NXF_CONTAINER_ENTRYPOINT_OVERRIDE: 'true']  | true

    }

    @Unroll
    def 'should validate oci auto-pull mode' () {

        expect:
        CONFIG.canRunOciImage() == OCI_AUTO_PULL

        where:
        CONFIG                                      | OCI_AUTO_PULL
        new SingularityConfig([:])                  | false
        new SingularityConfig(ociAutoPull:false)    | false
        and:
        new SingularityConfig(ociMode:true)         | true
        new ApptainerConfig(ociMode:true)           | false
        and:
        new SingularityConfig(ociAutoPull:true)     | true
        new ApptainerConfig(ociAutoPull:true)       | true

    }

    def 'should get fusion options' () {
        expect:
        CONFIG.getFusionOptions() == EXPECTED

        where:
        CONFIG                                          | EXPECTED
        new DockerConfig([:])                           | '--rm --privileged'
        new PodmanConfig([:])                           | '--rm --privileged'
        and:
        new SingularityConfig([:])                      | null
        new SingularityConfig(ociMode:true)             | '-B /dev/fuse'
        new SingularityConfig(ociAutoPull:true)         | null
        new ApptainerConfig(oci:true)                   | null
    }

}
