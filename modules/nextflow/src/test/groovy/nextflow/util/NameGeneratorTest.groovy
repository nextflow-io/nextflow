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

package nextflow.util

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NameGeneratorTest extends Specification {

    def 'should return a random name' () {
        when:
        def (adj, name)= NameGenerator.next().tokenize('_')
        then:
        NameGenerator.ADJECTIVES.contains(adj)
        NameGenerator.NAMES.contains(name)

    }

    def 'should not generate a random name except the specified one' () {
        when:
        def name = NameGenerator.next('evil_pike')
        then:
        name
        name != 'evil_pike'
    }
}
