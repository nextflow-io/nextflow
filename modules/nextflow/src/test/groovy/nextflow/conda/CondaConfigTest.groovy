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
 *
 */

package nextflow.conda

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondaConfigTest extends Specification {

    @Unroll
    def 'should check enabled flag'() {
        given:
        def conda = new CondaConfig(CONFIG, ENV)
        expect:
        conda.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]
        and:
        false       | [:]               | [NXF_CONDA_ENABLED: false]
        true        | [:]               | [NXF_CONDA_ENABLED: true]
        false       | [enabled: false]  | [NXF_CONDA_ENABLED: true]  // <-- config has priority
        true        | [enabled: true]   | [NXF_CONDA_ENABLED: true]
    }


    @Unroll
    def 'should check channels options'() {
        given:
        def conda = new CondaConfig(CONFIG, [:])
        expect:
        conda.getChannels() == EXPECTED

        where:
        EXPECTED                | CONFIG
        []                      | [:]
        ['bioconda']            | [channels:'bioconda']
        ['bioconda']            | [channels:['bioconda']]
        and:
        ['this','that','other'] | [channels:'this,that,other']
        ['this','that','other'] | [channels:'this, that ,other']
        and:
        ['this','that','other'] | [channels:['this','that','other']]
    }
}
