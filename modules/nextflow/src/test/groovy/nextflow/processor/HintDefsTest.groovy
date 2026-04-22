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

/**
 * Tests for {@link HintDefs}
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HintDefsTest extends Specification {

    def 'should accept valid hints'() {
        when:
        HintDefs.validateHints([
            consumableResources: 'my-license=1',
            'awsbatch/consumableResources': 'a=1,b=2',
            'seqera/machineRequirement.provisioning': 'spot',
        ])
        then:
        noExceptionThrown()
    }

    def 'should accept null and empty maps'() {
        expect:
        HintDefs.validateHints(null) == null
        HintDefs.validateHints([:]) == null
    }

    def 'should accept null hint value'() {
        when:
        HintDefs.validateHints([consumableResources: null])
        then:
        noExceptionThrown()
    }

    def 'should reject non-string value'() {
        when:
        HintDefs.validateHints([consumableResources: 42])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("Invalid hint value")
        e.message.contains("consumableResources")
    }

    def 'should reject list value'() {
        when:
        HintDefs.validateHints([consumableResources: ['a', 'b']])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject closure value'() {
        when:
        HintDefs.validateHints([consumableResources: { 'x' }])
        then:
        thrown(IllegalArgumentException)
    }

    def 'should reject empty key'() {
        when:
        HintDefs.validateHints(['': 'x'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("null or empty")
    }

    def 'should reject multi-segment key'() {
        when:
        HintDefs.validateHints(['a/b/c': 'x'])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains("a/b/c")
    }

}
