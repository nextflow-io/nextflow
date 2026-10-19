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

package nextflow.platform

import nextflow.NextflowMeta
import nextflow.config.Manifest
import nextflow.script.PlatformMetadata
import nextflow.script.WorkflowMetadata
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Tests for the auto resource labels helper
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AutoLabelsTest extends Specification {

    def 'should hold all the valid names in order'() {
        expect:
        AutoLabels.VALID_NAMES as List == [
            'projectName', 'userName', 'runName', 'sessionId', 'resume',
            'revision', 'commitId', 'repository', 'manifestName',
            'runtimeVersion', 'workflowId', 'workspaceId', 'computeEnvId'
        ]
    }

    def 'should parse a disabled value'() {
        expect:
        AutoLabels.parse(null) == [] as Set
        AutoLabels.parse(false) == [] as Set
    }

    def 'should parse a true value as all names'() {
        expect:
        AutoLabels.parse(true) == AutoLabels.VALID_NAMES
    }

    @Unroll
    def 'should parse a subset of names'() {
        expect:
        AutoLabels.parse(VALUE) == EXPECTED as Set

        where:
        VALUE                           | EXPECTED
        'runName'                       | ['runName']
        'runName,projectName'           | ['runName', 'projectName']
        ' runName , projectName '       | ['runName', 'projectName']
        'runName,,projectName'          | ['runName', 'projectName']
        ['runName', 'projectName']      | ['runName', 'projectName']
        [' runName ', null]             | ['runName']
        []                              | []
        ''                              | []
    }

    def 'should reject an unknown name'() {
        when:
        AutoLabels.parse('runName,foo')
        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('foo')
        e.message.contains('valid names are')
    }

    def 'should report the option name in the error message'() {
        when:
        AutoLabels.parse('foo', 'tower.autoLabels')
        then:
        def e1 = thrown(IllegalArgumentException)
        e1.message.contains("'tower.autoLabels'")

        when:
        AutoLabels.parse(new Object(), 'seqera.executor.autoLabels')
        then:
        def e2 = thrown(IllegalArgumentException)
        e2.message.contains("'seqera.executor.autoLabels'")
        e2.message.contains('expected true, false, a list, or a comma-separated string')
    }

    def 'should reject an unsupported value type'() {
        when:
        AutoLabels.parse(123)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should map all workflow metadata to labels'() {
        given:
        def sessionId = UUID.randomUUID()
        def platform = new PlatformMetadata('wf-abc123')
        platform.workspace = new PlatformMetadata.Workspace(workspaceId: '1234')
        platform.computeEnv = new PlatformMetadata.ComputeEnv(id: 'ce-abc')
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
            getPlatform() >> platform
        }

        when:
        def labels = AutoLabels.labelsFor(workflow, AutoLabels.VALID_NAMES)

        then:
        labels == [
            'nextflow.io/projectName': 'nf-core/rnaseq',
            'nextflow.io/userName': 'pditommaso',
            'nextflow.io/runName': 'crazy_darwin',
            'nextflow.io/sessionId': sessionId.toString(),
            'nextflow.io/resume': 'true',
            'nextflow.io/revision': '3.12.0',
            'nextflow.io/commitId': 'abc1234',
            'nextflow.io/repository': 'https://github.com/nf-core/rnaseq',
            'nextflow.io/manifestName': 'nf-core/rnaseq',
            'nextflow.io/runtimeVersion': NextflowMeta.instance.version.toString(),
            'seqera.io/platform/workflowId': 'wf-abc123',
            'seqera.io/platform/workspaceId': '1234',
            'seqera.io/platform/computeEnvId': 'ce-abc'
        ]
    }

    @Unroll
    def 'should prefer the platform user name over the OS user name'() {
        given:
        def platform = new PlatformMetadata()
        if( PLATFORM_USER )
            platform.user = new PlatformMetadata.User(id: '1', userName: PLATFORM_USER)
        def workflow = Mock(WorkflowMetadata) {
            getUserName() >> 'ec2-user'
            getPlatform() >> platform
        }

        expect:
        AutoLabels.labelsFor(workflow, ['userName'] as Set) == ['nextflow.io/userName': EXPECTED]

        where:
        PLATFORM_USER   | EXPECTED
        'jane.doe'      | 'jane.doe'
        null            | 'ec2-user'
    }

    def 'should omit the labels whose metadata is absent'() {
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
            getPlatform() >> new PlatformMetadata()
        }

        when:
        def labels = AutoLabels.labelsFor(workflow, AutoLabels.VALID_NAMES)

        then:
        labels.containsKey('nextflow.io/projectName')
        labels.containsKey('nextflow.io/userName')
        labels.containsKey('nextflow.io/runName')
        labels.containsKey('nextflow.io/sessionId')
        labels['nextflow.io/resume'] == 'false'
        and:
        !labels.containsKey('nextflow.io/revision')
        !labels.containsKey('nextflow.io/commitId')
        !labels.containsKey('nextflow.io/repository')
        !labels.containsKey('nextflow.io/manifestName')
        !labels.containsKey('seqera.io/platform/workflowId')
        !labels.containsKey('seqera.io/platform/workspaceId')
        !labels.containsKey('seqera.io/platform/computeEnvId')
        and:
        // no entry ever carries an empty value
        labels.values().every { it }
    }

    def 'should omit the platform labels when the metadata is missing altogether'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getRunName() >> 'happy_turing'
            isResume() >> false
            getManifest() >> new Manifest([:])
            getPlatform() >> null
        }

        when:
        def labels = AutoLabels.labelsFor(workflow, AutoLabels.VALID_NAMES)

        then:
        labels.keySet() == ['nextflow.io/runName', 'nextflow.io/resume', 'nextflow.io/runtimeVersion'] as Set
    }

    def 'should map only the included names'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'nf-core/rnaseq'
            getRunName() >> 'crazy_darwin'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getRevision() >> '3.12.0'
            getManifest() >> new Manifest([name: 'nf-core/rnaseq'])
        }

        expect:
        AutoLabels.labelsFor(workflow, ['runName', 'revision'] as Set).keySet() == ['nextflow.io/runName', 'nextflow.io/revision'] as Set
    }

    def 'should map nothing for an empty include set or missing metadata'() {
        given:
        def workflow = Mock(WorkflowMetadata) { getProjectName() >> 'hello' }

        expect:
        AutoLabels.labelsFor(workflow, [] as Set) == [:]
        AutoLabels.labelsFor(workflow, null) == [:]
        AutoLabels.labelsFor(null, AutoLabels.VALID_NAMES) == [:]
    }

}
