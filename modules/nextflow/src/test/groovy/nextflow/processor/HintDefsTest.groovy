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

    def 'should accept known hint key'() {
        when:
        HintDefs.validateHints([consumableResources: 'my-license'])
        then:
        noExceptionThrown()
    }

    def 'should warn on unknown key with close match'() {
        // consumableResource is close to consumableResources
        when:
        HintDefs.validateHints([consumableResource: 'my-license'])
        then:
        noExceptionThrown()
        // warning is logged — verified by log output in integration tests
    }

    def 'should warn on unknown key with no close match'() {
        when:
        HintDefs.validateHints([somethingRandom: 'value'])
        then:
        noExceptionThrown()
    }

    def 'should reject invalid value type'() {
        when:
        HintDefs.validateHints([consumableResources: ['a', 'b']])
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Invalid hint value type')
        e.message.contains('consumableResources')
    }

    def 'should accept string and integer values'() {
        when:
        HintDefs.validateHints([consumableResources: 'my-license', 'scheduling.priority': 10])
        then:
        noExceptionThrown()
    }

    def 'should skip executor-prefixed keys'() {
        when:
        HintDefs.validateHints(['seqera/machineRequirement.arch': 'arm64', 'seqera/unknownKey': 'value'])
        then:
        noExceptionThrown()
    }

    def 'should handle null and empty maps'() {
        when:
        HintDefs.validateHints(null)
        then:
        noExceptionThrown()

        when:
        HintDefs.validateHints([:])
        then:
        noExceptionThrown()
    }

    def 'should accept null hint value'() {
        when:
        HintDefs.validateHints([consumableResources: null])
        then:
        noExceptionThrown()
    }

}
