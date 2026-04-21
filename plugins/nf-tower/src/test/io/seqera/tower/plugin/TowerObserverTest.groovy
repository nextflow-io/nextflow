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

package io.seqera.tower.plugin

import nextflow.script.PlatformMetadata

import java.net.http.HttpResponse
import java.nio.file.Files
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId

import com.github.tomakehurst.wiremock.WireMockServer
import com.github.tomakehurst.wiremock.client.WireMock
import io.seqera.http.HxClient
import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.container.DockerConfig
import nextflow.container.resolver.ContainerMeta
import nextflow.exception.AbortOperationException
import nextflow.script.PlatformMetadata
import nextflow.script.ScriptBinding
import nextflow.script.WorkflowMetadata
import nextflow.trace.TraceRecord
import nextflow.trace.WorkflowStats
import nextflow.trace.WorkflowStatsObserver
import nextflow.util.Duration
import nextflow.util.ProcessHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerObserverTest extends Specification {

    protected boolean aroundNow(value) {
        def now = Instant.now().toEpochMilli()
        value > now-1_000 && value <= now
    }

    private TowerObserver newObserver(Session session, Map env = [:]) {
        def client = Mock(TowerClient)
        def observer = new TowerObserver(session, client, null, env)
        observer.@reports = Mock(TowerReports)
        return observer
    }

    def 'should create message map' () {
        given:
        def session = Mock(Session)
        def params = new ScriptBinding.ParamsMap(x: "hello")
        def meta = Mock(WorkflowMetadata)
        def observer = Spy(newObserver(session))
        observer.@workflowId = '12ef'

        when:
        def map = observer.makeCompleteReq(session)
        then:
        1 * session.getWorkflowMetadata() >> meta
        1 * session.getParams() >> params
        1 * meta.toMap() >> [foo:1, bar:2, container: [p1: 'c1', p2: 'c2']]
        1 * observer.getMetricsList() >> [[process:'foo', cpu: [min: 1, max:5], time: [min: 6, max: 9]]]
        1 * observer.getWorkflowProgress(false) >> new WorkflowProgress()
        1 * observer.getOutFile() >> 'bar.out'
        1 * observer.getLogFile() >> 'foo.out'
        1 * observer.getOperationId() >> 'op-12345'
        then:
        map.workflow.foo == 1
        map.workflow.bar == 2
        map.workflow.id == '12ef'
        map.workflow.params == [x: 'hello']
        map.workflow.container == null
        map.metrics == [[process:'foo', cpu: [min: 1, max:5], time: [min: 6, max: 9]]]
        map.progress == new WorkflowProgress()
        and:
        aroundNow(map.instant)
        and:
        map.workflow.outFile == 'bar.out'
        map.workflow.logFile == 'foo.out'
        map.workflow.operationId == 'op-12345'
    }

    def 'should capitalise underscores' () {
        given:
        def tower = new TowerObserver(Mock(Session), Mock(TowerClient), "ws1234", [:] )

        expect:
        tower.underscoreToCamelCase(STR) == EXPECTED
        where:
        STR         | EXPECTED
        'abc'       | 'abc'
        'a_b_c'     | 'aBC'
        'foo__bar'  | 'fooBar'
    }

    def 'should post task records' () {
        given:
        def session = Mock(Session)
        def PROGRESS = Mock(WorkflowProgress) { getRunning()>>1; getSucceeded()>>2; getFailed()>>3 }
        def observer = Spy(newObserver(session))
        observer.@workflowId = 'xyz-123'

        def nowTs = System.currentTimeMillis()
        def submitTs = nowTs-2000
        def startTs = nowTs-1000

        def trace = new TraceRecord([
                taskId: 10,
                process: 'foo',
                workdir: "/work/dir",
                cpus: 1,
                submit: submitTs,
                start: startTs,
                complete: nowTs ])
        trace.executorName= 'batch'
        trace.machineInfo = new CloudMachineInfo('m4.large', 'eu-west-1b', PriceModel.spot)
        trace.containerMeta = new ContainerMeta(requestId: '12345', sourceImage: 'ubuntu:latest', targetImage: 'wave.io/12345/ubuntu:latest')

        when:
        def req = observer.makeTasksReq([trace])
        then:
        observer.getWorkflowProgress(true) >> PROGRESS
        and:
        req.tasks[0].taskId == 10
        req.tasks[0].process == 'foo'
        req.tasks[0].workdir == "/work/dir"
        req.tasks[0].cpus == 1
        req.tasks[0].submit == OffsetDateTime.ofInstant(Instant.ofEpochMilli(submitTs), ZoneId.systemDefault())
        req.tasks[0].start == OffsetDateTime.ofInstant(Instant.ofEpochMilli(startTs), ZoneId.systemDefault())
        req.tasks[0].executor == 'batch'
        req.tasks[0].machineType == 'm4.large'
        req.tasks[0].cloudZone == 'eu-west-1b'
        req.tasks[0].priceModel == 'spot'
        and:
        req.progress.running == 1
        req.progress.succeeded == 2
        req.progress.failed == 3
        and:
        req.containers[0].requestId == '12345'
        req.containers[0].sourceImage == 'ubuntu:latest'
        req.containers[0].targetImage == 'wave.io/12345/ubuntu:latest'
        and:
        aroundNow(req.instant)
    }

    static now_millis = System.currentTimeMillis()
    static now_instant = OffsetDateTime.ofInstant(Instant.ofEpochMilli(now_millis), ZoneId.systemDefault())

    def 'should fix field types' () {

        expect:
        TowerObserver.fixTaskField(FIELD, VALUE) == EXPECTED

        where:
        FIELD       | VALUE         | EXPECTED
        'foo'       | 'hola'        | 'hola'
        'submit'    | now_millis    | now_instant
        'start'     | now_millis    | now_instant
        'complete'  | now_millis    | now_instant
        'complete'  | 0             | null
    }

    def 'should create workflow json' () {

        given:
        def sessionId = UUID.randomUUID()
        def dir = Files.createTempDirectory('test')
        def session = Mock(Session)
        session.getUniqueId() >> sessionId
        session.getRunName() >> 'foo'
        session.config >> [:]
        session.containerConfig >> new DockerConfig([:])
        session.getParams() >> new ScriptBinding.ParamsMap([foo:'Hello', bar:'World'])

        def meta = new WorkflowMetadata(
                session: session,
                projectName: 'the-project-name',
                repository: 'git://repo.com/foo' )
        session.getWorkflowMetadata() >> meta
        session.getStatsObserver() >> Mock(WorkflowStatsObserver) { getStats() >> new WorkflowStats() }

        def observer = Spy(newObserver(session, ENV))
        observer.getOperationId() >> 'op-112233'
        observer.getLogFile() >> 'log.file'
        observer.getOutFile() >> 'out.file'

        when:
        def req1 = observer.makeCreateReq(session)
        then:
        req1.sessionId == sessionId.toString()
        req1.runName == 'foo'
        req1.projectName == 'the-project-name'
        req1.repository == 'git://repo.com/foo'
        req1.workflowId == WORKFLOW_ID
        and:
        aroundNow(req1.instant)

        when:
        def req = observer.makeBeginReq(session)
        then:
        observer.getWorkflowId() >> '12345'
        and:
        req.workflow.id == '12345'
        req.workflow.params == [foo:'Hello', bar:'World']
        req.workflow.outFile == 'out.file'
        req.workflow.logFile == 'log.file'
        req.workflow.operationId == 'op-112233'
        and:
        req.towerLaunch == TOWER_LAUNCH
        and:
        aroundNow(req.instant)

        cleanup:
        dir?.deleteDir()

        where:
        ENV                         | WORKFLOW_ID   | TOWER_LAUNCH
        [:]                         | null          | false
        [TOWER_WORKFLOW_ID: '1234'] | '1234'        | true

    }

    def 'should convert map' () {
        given:
        def tower = new TowerObserver(Mock(Session), Mock(TowerClient), "ws1234", [:] )

        expect:
        tower.mapToString(null)  == null
        tower.mapToString('ciao') == 'ciao'
        tower.mapToString([:]) == null
        tower.mapToString([p:'foo', q:'bar']) == null
    }

    def 'should create init request' () {
        given:
        def uuid = UUID.randomUUID()
        def meta = Mock(WorkflowMetadata) {
            getProjectName() >> 'the-project-name'
            getRepository() >> 'git://repo.com/foo'
        }
        def session = Mock(Session) {
            getUniqueId() >> uuid
            getRunName() >> 'foo_bar'
            getWorkflowMetadata() >> meta
        }
        def observer = newObserver(session, [TOWER_WORKFLOW_ID: 'x123'])

        when:
        def req = observer.makeCreateReq(session)
        then:
        req.sessionId == uuid.toString()
        req.runName == 'foo_bar'
        req.projectName == 'the-project-name'
        req.repository == 'git://repo.com/foo'
        req.workflowId == 'x123'
        and:
        aroundNow(req.instant)

        and:
        observer.towerLaunch
    }

    def 'should post create request' () {
        given:
        def uuid = UUID.randomUUID()
        def platform = new PlatformMetadata()
        def meta = Mock(WorkflowMetadata) {
            getProjectName() >> 'the-project-name'
            getRepository() >> 'git://repo.com/foo'
            getPlatform() >> platform
        }
        def session = Mock(Session) {
            getUniqueId() >> uuid
            getRunName() >> 'foo_bar'
            getWorkflowMetadata() >> meta
        }
        def towerClient = Mock(TowerClient)
        def observer = Spy(new TowerObserver(session, towerClient, null, [:]))
        observer.@reports = Mock(TowerReports)

        when:
        observer.onFlowCreate(session)
        then:
        1 * observer.makeCreateReq(session) >> [runName: 'foo']
        1 * towerClient.traceCreate([runName: 'foo'], null) >> [workflowId: 'xyz123', watchUrl: 'https://cloud.seqera.io/watch/xyz123']
        and:
        observer.runName == 'foo_bar'
        observer.runId == uuid.toString()
        and:
        observer.workflowId == 'xyz123'
        observer.@watchUrl == 'https://cloud.seqera.io/watch/xyz123'
        !observer.towerLaunch
        and:
        platform.workflowId == 'xyz123'
        platform.workflowUrl == 'https://cloud.seqera.io/watch/xyz123'

    }

    def 'should set workflowUrl on platform metadata during onFlowBegin' () {
        given:
        def platform = new PlatformMetadata()
        def meta = Mock(WorkflowMetadata) {
            getPlatform() >> platform
        }
        def session = Mock(Session) {
            getWorkflowMetadata() >> meta
        }
        def towerClient = Mock(TowerClient)
        def observer = Spy(new TowerObserver(session, towerClient, null, [:]))
        observer.@reports = Mock(TowerReports)
        observer.@workflowId = 'abc123'

        when:
        observer.onFlowBegin()
        then:
        1 * observer.makeBeginReq(session) >> [foo: 'bar']
        1 * towerClient.traceBegin([foo: 'bar'], null, 'abc123') >> [watchUrl: 'https://cloud.seqera.io/watch/abc123']
        and:
        observer.@watchUrl == 'https://cloud.seqera.io/watch/abc123'
        platform.workflowUrl == 'https://cloud.seqera.io/watch/abc123'

        cleanup:
        observer.@sender?.interrupt()
    }

    def 'should fetch workflow meta' () {
        given:
        def session = Mock(Session)
        def observer = newObserver(session, ENV)

        expect:
        observer.getOperationId() == OP_ID
        observer.getLogFile() == LOG_FILE
        observer.getOutFile() == OUT_FILE

        where:
        OP_ID                                           | OUT_FILE      | LOG_FILE    | ENV
        null                                            | null          | null        | [:]
        "local-platform::${ProcessHelper.selfPid()}"    | null          | null        | [TOWER_ALLOW_NEXTFLOW_LOGS:'true']
        'aws-batch::1234z'                              | 'xyz.out'     | 'hola.log'  | [TOWER_ALLOW_NEXTFLOW_LOGS:'true', AWS_BATCH_JOB_ID: '1234z', NXF_OUT_FILE: 'xyz.out', NXF_LOG_FILE: 'hola.log']
    }

    def 'should deduplicate containers' () {
        given:
        def session = Mock(Session)
        def observer = newObserver(session)
        and:
        def c1 = new ContainerMeta(requestId: '12345', sourceImage: 'ubuntu:latest', targetImage: 'wave.io/12345/ubuntu:latest')
        def c2 = new ContainerMeta(requestId: '54321', sourceImage: 'ubuntu:latest', targetImage: 'wave.io/54321/ubuntu:latest')
        and:
        def trace1 = new TraceRecord(
                taskId: 1,
                process: 'foo',
                workdir: "/work/dir",
                cpus: 1,
                submit: System.currentTimeMillis(),
                start: System.currentTimeMillis(),
                complete: System.currentTimeMillis())
        trace1.containerMeta = c1
        and:
        def trace2 = new TraceRecord(
            taskId: 2,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: System.currentTimeMillis(),
            start: System.currentTimeMillis(),
            complete: System.currentTimeMillis())
        trace2.containerMeta = c2
        and:
        def trace3 = new TraceRecord(
            taskId: 3,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: System.currentTimeMillis(),
            start: System.currentTimeMillis(),
            complete: System.currentTimeMillis())
        trace3.containerMeta = c2

        expect:
        observer.getNewContainers([trace1]) == [c1]
        and:
        observer.getNewContainers([trace1]) == []
        and:
        observer.getNewContainers([trace1, trace2, trace3]) == [c2]
    }

    def 'should not send complete request when onFlowBegin was not invoked' () {
        given:
        def session = Mock(Session)
        def towerClient = Mock(TowerClient)
        def observer = Spy(new TowerObserver(session, towerClient, null, [:]))
        def reports = Mock(TowerReports)
        observer.@reports = reports
        observer.@workflowId = 'xyz-123'
        observer.@sender = null

        when:
        observer.onFlowComplete()

        then:
        1 * reports.publishRuntimeReports()
        1 * reports.flowComplete()
        0 * towerClient.traceComplete(_, _, _)
    }

    def 'should apply platform metadata from trace create response'() {
        given:
        def metadata = new WorkflowMetadata()
        def session = Mock(Session) {
            getWorkflowMetadata() >> metadata
        }
        def observer = new TowerObserver(session, Mock(TowerClient), '1234',  SysEnv.get())

        def responseMetadata = [
            userId: 39,
            userName: 'user',
            userOrganization: 'ACME Inc.',
            workspaceId: 1234,
            workspaceName: 'Workspace-Name',
            workspaceFullName: 'Full Workspace Name',
            orgName: 'ACME Inc.',
            computeEnvId: 'ce1234',
            computeEnvName: 'ce-test',
            computeEnvPlatform: 'aws-batch',
            pipelineName: 'test-pipeline',
            pipelineId: 'pipe1234',
            revision: 'v1.1',
            commitId: 'abcd12345'
        ]

        when:
        observer.applyPlatformMetadata(responseMetadata)

        then:
        metadata.platform.user.id == '39'
        metadata.platform.user.userName == 'user'
        metadata.platform.user.organization == 'ACME Inc.'
        metadata.platform.workspace.id == '1234'
        metadata.platform.workspace.name == 'Workspace-Name'
        metadata.platform.workspace.fullName == 'Full Workspace Name'
        metadata.platform.workspace.organization == 'ACME Inc.'
        metadata.platform.computeEnv.id == 'ce1234'
        metadata.platform.computeEnv.name == 'ce-test'
        metadata.platform.computeEnv.platform == 'aws-batch'
        metadata.platform.pipeline.id == 'pipe1234'
        metadata.platform.pipeline.name == 'test-pipeline'
        metadata.platform.pipeline.revision == 'v1.1'
        metadata.platform.pipeline.commitId == 'abcd12345'
    }

    def 'should include numSpotInterruptions in task map'() {
        given:
        def session = Mock(Session)
        def observer = Spy(newObserver(session))
        observer.getWorkflowProgress(true) >> new WorkflowProgress()

        def now = System.currentTimeMillis()
        def trace = new TraceRecord([
            taskId: 42,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: now-2000,
            start: now-1000,
            complete: now
        ])
        trace.setNumSpotInterruptions(3)

        when:
        def req = observer.makeTasksReq([trace])

        then:
        req.tasks.size() == 1
        req.tasks[0].numSpotInterruptions == 3
    }

    def 'should include logStreamId in task map'() {
        given:
        def session = Mock(Session)
        def observer = Spy(newObserver(session))
        observer.getWorkflowProgress(true) >> new WorkflowProgress()

        def now = System.currentTimeMillis()
        def trace = new TraceRecord([
            taskId: 42,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: now-2000,
            start: now-1000,
            complete: now
        ])
        trace.setLogStreamId('arn:aws:logs:us-east-1:123456789:log-group:/ecs/task:log-stream:abc123')

        when:
        def req = observer.makeTasksReq([trace])

        then:
        req.tasks.size() == 1
        req.tasks[0].logStreamId == 'arn:aws:logs:us-east-1:123456789:log-group:/ecs/task:log-stream:abc123'
    }

    def 'should include resourceAllocation in task map'() {
        given:
        def session = Mock(Session)
        def observer = Spy(newObserver(session))
        observer.getWorkflowProgress(true) >> new WorkflowProgress()

        def now = System.currentTimeMillis()
        def trace = new TraceRecord([
            taskId: 42,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: now-2000,
            start: now-1000,
            complete: now,
            accelerator: 2,
            acceleratorType: 'v100'
        ])
        trace.setResourceAllocation([cpuShares: 2048, memoryMiB: 4096, time: '1h'])

        when:
        def req = observer.makeTasksReq([trace])

        then:
        req.tasks.size() == 1
        req.tasks[0].accelerator == 2
        req.tasks[0].acceleratorType == 'v100'
        req.tasks[0].resourceAllocation == [cpuShares: 2048, memoryMiB: 4096, time: '1h']
    }

    def 'should include gpuMetrics in task map'() {
        given:
        def session = Mock(Session)
        def observer = Spy(newObserver(session))
        observer.getWorkflowProgress(true) >> new WorkflowProgress()

        def now = System.currentTimeMillis()
        def trace = new TraceRecord([
            taskId: 42,
            process: 'foo',
            workdir: "/work/dir",
            cpus: 1,
            submit: now-2000,
            start: now-1000,
            complete: now
        ])
        trace.setGpuMetrics([name: 'Tesla T4', mem: 15360, driver: '580.126.09', active_time: 651030, pct: 75, peak: 100])

        when:
        def req = observer.makeTasksReq([trace])

        then:
        req.tasks.size() == 1
        req.tasks[0].gpuMetrics.name == 'Tesla T4'
        req.tasks[0].gpuMetrics.mem == 15360
        req.tasks[0].gpuMetrics.pct == 75
        req.tasks[0].gpuMetrics.peak == 100
    }


}
