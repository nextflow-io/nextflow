/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

}
