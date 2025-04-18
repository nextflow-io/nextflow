/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage

import nextflow.lineage.config.LineageConfig
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DefaultLinStoreFactoryTest extends Specification {

    @Unroll
    def 'should validate can open' () {
        given:
        def factory = new DefaultLinStoreFactory()
        def config = new LineageConfig(CONFIG)

        expect:
        factory.canOpen(config) == EXPECTED
        
        where:
        EXPECTED    | CONFIG
        true        | [:]
        true        | [store:[location:'/some/path']]
        true        | [store:[location:'some/rel/path']]
        true        | [store:[location:'file:/this/that']]
        true        | [store:[location:'s3://some/path']]
        false       | [store:[location:'http://some/path']]
        false       | [store:[location:'jdbc:foo']]
    }

}
