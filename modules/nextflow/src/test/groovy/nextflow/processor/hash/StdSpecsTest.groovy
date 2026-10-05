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
package nextflow.processor.hash

import nextflow.util.EncodingRules
import spock.lang.Specification

class StdSpecsTest extends Specification {

    def 'adjacent specs differ in exactly the documented place'() {
        expect: 'v1 to v2 differs only in encoding'
        StdSpecs.STD_V1_3.keys() == StdSpecs.STD_V1_5.keys()
        StdSpecs.STD_V1_3.encoding.canonicalForm() == EncodingRules.LEGACY.canonicalForm()
        StdSpecs.STD_V1_5.encoding.canonicalForm() == EncodingRules.RECORD_TYPES.canonicalForm()

        and: 'v2 to v3 differs only by the module bundle key'
        StdSpecs.STD_V1_6.keys() - StdSpecs.STD_V1_5.keys() == [HashKey.RESOURCES_BUNDLE]
        StdSpecs.STD_V1_5.keys() - StdSpecs.STD_V1_6.keys() == []
        StdSpecs.STD_V1_6.encoding.canonicalForm() == StdSpecs.STD_V1_5.encoding.canonicalForm()

        and: 'v3 to v4 differs only in the eval extraction'
        StdSpecs.STD_V1_6.keys() == StdSpecs.STD_V1_7.keys()
        StdSpecs.STD_V1_6.bindings[HashKey.EVAL_OUTPUTS].canonicalName() ==
            'evalOutputs.derivedString'
        StdSpecs.STD_V1_7.bindings[HashKey.EVAL_OUTPUTS].canonicalName() ==
            'evalOutputs.rawMap'
    }

    def 'bin entries are present in every spec, as they have been since 2018'() {
        expect:
        StdSpecs.all().every { HashKey.BIN_ENTRIES in it.keys() }
    }

    def 'every spec accounts for every key in the vocabulary'() {
        given: 'keys a spec may legitimately omit, with the reason'
        def permittedOmissions = [
            'std/v1.1': [HashKey.RESOURCES_BUNDLE] as Set,
            'std/v1.2': [HashKey.RESOURCES_BUNDLE] as Set,
            'std/v1.3': [HashKey.RESOURCES_BUNDLE] as Set,   // #6914 (2026-07-17) ends v1.5
            'std/v1.4': [HashKey.RESOURCES_BUNDLE] as Set,
            'std/v1.5': [HashKey.RESOURCES_BUNDLE] as Set,
            'std/v1.6': [] as Set,
            'std/v1.7': [] as Set
        ]

        expect:
        StdSpecs.all().every { spec ->
            def missing = (HashKey.values() as Set) - (spec.keys() as Set)
            missing == permittedOmissions[spec.id]
        }
    }

    def 'fingerprints are distinct across all specs'() {
        expect:
        StdSpecs.all()*.fingerprint().toSet().size() == StdSpecs.all().size()
    }
}
