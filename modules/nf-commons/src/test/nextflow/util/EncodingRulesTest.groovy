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

package nextflow.util

import com.google.common.hash.Hashing
import spock.lang.Specification

class EncodingRulesTest extends Specification {

    def 'legacy rules hash a map by values, record-types rules by entries'() {
        given:
        def map = [alpha: 'one', beta: 'two']

        when:
        def legacy = EncodingRules.LEGACY
            .apply(new HashBuilder().withHasher(Hashing.murmur3_128().newHasher()))
            .with(map)
            .build()
        and:
        def recordTypes = EncodingRules.RECORD_TYPES
            .apply(new HashBuilder().withHasher(Hashing.murmur3_128().newHasher()))
            .with(map)
            .build()

        then:
        legacy != recordTypes
    }

    def 'legacy rules hash only the values, so renaming a key does not change the hash'() {
        given:
        def one = [alpha: 'x']
        def two = [renamed: 'x']

        expect:
        EncodingRules.LEGACY.apply(new HashBuilder()).with(one).build() ==
            EncodingRules.LEGACY.apply(new HashBuilder()).with(two).build()
        and:
        EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(one).build() !=
            EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(two).build()
    }

    def 'default HashBuilder behaviour is unchanged'() {
        given:
        def map = [alpha: 'one', beta: 'two']

        expect:
        new HashBuilder().with(map).build() ==
            EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(map).build()
    }

    def 'canonical form is stable and distinguishes the rule sets'() {
        expect:
        EncodingRules.LEGACY.canonicalForm() == 'orderIndependentMaps=false;cacheFunnelFirst=false'
        and:
        EncodingRules.RECORD_TYPES.canonicalForm() == 'orderIndependentMaps=true;cacheFunnelFirst=true'
    }
}
