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

package io.seqera.executor

import spock.lang.Specification

/**
 * Tests for Labels helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LabelsTest extends Specification {

    def 'should merge label maps coercing values with the overlay winning on collision'() {
        expect:
        Labels.merge(
                ['nextflow.io/runName': 'happy_turing', 'nextflow.io/projectName': 'hello'],
                ['nextflow.io/runName': 'custom', team: 7]
        ) == ['nextflow.io/runName': 'custom', 'nextflow.io/projectName': 'hello', team: '7']
    }

    def 'should treat a null or empty map as no labels when merging'() {
        expect:
        Labels.merge(null, null) == [:]
        Labels.merge(['team': 'a'], null) == ['team': 'a']
        Labels.merge(null, ['team': 'a']) == ['team': 'a']
    }

    def 'should coerce map values to strings'() {
        expect:
        Labels.toStringMap(null) == [:]
        Labels.toStringMap([:]) == [:]
        Labels.toStringMap([a: 1, b: 'x', c: true]) == [a: '1', b: 'x', c: 'true']
    }

    def 'should reject non-map resourceLabels with a clear error'() {
        when:
        Labels.toStringMap(['foo', 'bar'])

        then:
        def err = thrown(IllegalArgumentException)
        err.message.contains("'resourceLabels'")
        err.message.contains('map of key/value pairs')
        err.message.contains('java.util.ArrayList')
    }

    def 'should compute null delta when task labels are empty'() {
        expect:
        Labels.delta(null, [team: 'a']) == null
        Labels.delta([:], [team: 'a']) == null
    }

    def 'should return full task labels when run labels are empty'() {
        expect:
        Labels.delta([team: 'a', region: 'us'], null) == [team: 'a', region: 'us']
        Labels.delta([team: 'a', region: 'us'], [:]) == [team: 'a', region: 'us']
    }

    def 'should keep only differing or missing keys in delta'() {
        expect:
        Labels.delta([team: 'a', region: 'us'], [team: 'a']) == [region: 'us']
        Labels.delta([team: 'b'], [team: 'a']) == [team: 'b']
        Labels.delta([team: 'a', region: 'us'], [team: 'a', region: 'us']) == null
    }
}
