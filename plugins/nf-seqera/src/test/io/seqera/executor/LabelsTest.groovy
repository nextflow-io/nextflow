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

import nextflow.NextflowMeta
import nextflow.config.Manifest
import nextflow.script.PlatformMetadata
import nextflow.script.WorkflowMetadata
import spock.lang.Specification

/**
 * Tests for Labels helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LabelsTest extends Specification {

    def 'should create labels with all workflow metadata'() {
        given:
        def sessionId = UUID.randomUUID()
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'nf-core/rnaseq'
            getUserName() >> 'pditommaso'
            getRunName() >> 'crazy_darwin'
            getSessionId() >> sessionId
            isResume() >> true
            getRevision() >> '3.12.0'
            getCommitId() >> 'abc1234'
            getRepository() >> 'https://github.com/nf-core/rnaseq'
            getManifest() >> new Manifest([name: 'nf-core/rnaseq'])
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)

        then:
        labels.entries['nextflow.io/projectName'] == 'nf-core/rnaseq'
        labels.entries['nextflow.io/userName'] == 'pditommaso'
        labels.entries['nextflow.io/runName'] == 'crazy_darwin'
        labels.entries['nextflow.io/sessionId'] == sessionId.toString()
        labels.entries['nextflow.io/resume'] == 'true'
        labels.entries['nextflow.io/revision'] == '3.12.0'
        labels.entries['nextflow.io/commitId'] == 'abc1234'
        labels.entries['nextflow.io/repository'] == 'https://github.com/nf-core/rnaseq'
        labels.entries['nextflow.io/manifestName'] == 'nf-core/rnaseq'
        labels.entries['nextflow.io/runtimeVersion'] == NextflowMeta.instance.version.toString()
    }

    def 'should compute stable runId from sessionId and runName'() {
        given:
        def sid = 'e2315a82-49b0-4langc3-a58a-0d7d52f7e3a1'
        def runName = 'crazy_darwin'

        expect:
        Labels.runId(sid, runName) == Labels.runId(sid, runName)
        Labels.runId(sid, runName) != Labels.runId(sid, 'other_name')
        Labels.runId(sid, runName) != Labels.runId(UUID.randomUUID().toString(), runName)
    }

    def 'should omit null workflow metadata from labels'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getUserName() >> 'user1'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getRevision() >> null
            getCommitId() >> null
            getRepository() >> null
            getManifest() >> new Manifest([:])
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)

        then:
        labels.entries.containsKey('nextflow.io/projectName')
        labels.entries.containsKey('nextflow.io/userName')
        labels.entries.containsKey('nextflow.io/runName')
        labels.entries.containsKey('nextflow.io/sessionId')
        labels.entries['nextflow.io/resume'] == 'false'
        !labels.entries.containsKey('nextflow.io/revision')
        !labels.entries.containsKey('nextflow.io/commitId')
        !labels.entries.containsKey('nextflow.io/repository')
        !labels.entries.containsKey('nextflow.io/manifestName')
    }

    def 'should add scheduler labels'() {
        when:
        def labels = new Labels()
                .withSchedRunId('run-123')
                .withSchedClusterId('cluster-456')

        then:
        labels.entries['seqera:sched:runId'] == 'run-123'
        labels.entries['seqera:sched:clusterId'] == 'cluster-456'
    }

    def 'should skip null scheduler labels'() {
        when:
        def labels = new Labels()
                .withSchedRunId(null)
                .withSchedClusterId(null)

        then:
        !labels.entries.containsKey('seqera:sched:runId')
        !labels.entries.containsKey('seqera:sched:clusterId')
    }

    def 'should allow user labels to override implicit labels'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getUserName() >> 'user1'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getManifest() >> new Manifest([:])
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)
                .withUserLabels([
                    'nextflow.io/runName': 'custom_name',
                    'team': 'research'
                ])

        then:
        labels.entries['nextflow.io/runName'] == 'custom_name'
        labels.entries['team'] == 'research'
        labels.entries['nextflow.io/projectName'] == 'hello'
    }

    def 'should include platform workflowId when available'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getUserName() >> 'user1'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getManifest() >> new Manifest([:])
            getPlatform() >> new PlatformMetadata('wf-abc123')
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)

        then:
        labels.entries['seqera.io/platform/workflowId'] == 'wf-abc123'
    }

    def 'should omit platform workflowId when not set'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getUserName() >> 'user1'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getManifest() >> new Manifest([:])
            getPlatform() >> new PlatformMetadata()
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)

        then:
        !labels.entries.containsKey('seqera.io/platform/workflowId')
    }

    def 'should handle null user labels'() {
        when:
        def labels = new Labels()
                .withUserLabels(null)

        then:
        labels.entries.isEmpty()
    }

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

    def 'should let process resource labels override workflow metadata on key collision'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getManifest() >> new Manifest([:])
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)
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
