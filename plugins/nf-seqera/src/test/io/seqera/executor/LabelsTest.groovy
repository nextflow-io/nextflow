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

    def 'should add process resource labels coercing values to string'() {
        when:
        def labels = new Labels()
                .withProcessResourceLabels([team: 'genomics', priority: 7, retain: true])

        then:
        labels.entries['team'] == 'genomics'
        labels.entries['priority'] == '7'
        labels.entries['retain'] == 'true'
    }

    def 'should ignore null or empty process resource labels'() {
        when:
        def a = new Labels().withProcessResourceLabels(null)
        def b = new Labels().withProcessResourceLabels([:])

        then:
        a.entries.isEmpty()
        b.entries.isEmpty()
    }

    def 'should let the resource labels added last win on key collision'() {
        when:
        def labels = new Labels()
                .withProcessResourceLabels(['nextflow.io/runName': 'happy_turing', 'nextflow.io/projectName': 'hello'])
                .withProcessResourceLabels(['nextflow.io/runName': 'custom', team: 'a'])

        then:
        labels.entries['nextflow.io/runName'] == 'custom'
        labels.entries['team'] == 'a'
        labels.entries['nextflow.io/projectName'] == 'hello'
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
