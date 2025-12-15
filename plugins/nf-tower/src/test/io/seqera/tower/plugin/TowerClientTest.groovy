/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package io.seqera.tower.plugin

import java.net.http.HttpResponse
import java.nio.file.Files
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId

import io.seqera.http.HxClient
import nextflow.Session
import nextflow.cloud.types.CloudMachineInfo
import nextflow.cloud.types.PriceModel
import nextflow.container.DockerConfig
import nextflow.container.resolver.ContainerMeta
import nextflow.exception.AbortOperationException
import nextflow.script.ScriptBinding
import nextflow.script.WorkflowMetadata
import nextflow.trace.TraceRecord
import nextflow.trace.WorkflowStats
import nextflow.trace.WorkflowStatsObserver
import nextflow.util.ProcessHelper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerClientTest extends Specification {

    protected boolean aroundNow(value) {
        def now = Instant.now().toEpochMilli()
        value > now-1_000 && value <= now
    }

    def 'should parse response' () {
        given:
        def tower = new TowerClient()

        when:
        def resp = new TowerClient.Response(200, '{"status":"OK", "workflowId":"12345", "watchUrl": "http://foo.com/watch/12345"}')
        def result = tower.parseTowerResponse(resp)
        then:
        result.workflowId == '12345'
        result.watchUrl == 'http://foo.com/watch/12345'

        when:
        resp = new TowerClient.Response(500, '{"status":"OK", "workflowId":"12345"}')
        tower.parseTowerResponse(resp)
        then:
        thrown(Exception)
    }

    def 'should create message map' () {
        given:
        def session = Mock(Session)
        def params = new ScriptBinding.ParamsMap(x: "hello")
        def meta = Mock(WorkflowMetadata)

        def tower = Spy(TowerClient)
        tower.@runName = session.runName
        tower.@workflowId = '12ef'

        when:
        def map = tower.makeCompleteReq(session)
        then:
        1 * session.getWorkflowMetadata() >> meta
        1 * session.getParams() >> params
        1 * meta.toMap() >> [foo:1, bar:2, container: [p1: 'c1', p2: 'c2']]
        1 * tower.getMetricsList() >> [[process:'foo', cpu: [min: 1, max:5], time: [min: 6, max: 9]]]
        1 * tower.getWorkflowProgress(false) >> new WorkflowProgress()
        1 * tower.getOutFile() >> 'bar.out'
        1 * tower.getLogFile() >> 'foo.out'
        1 * tower.getOperationId() >> 'op-12345'
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
        def tower = new TowerClient()

        expect:
        tower.underscoreToCamelCase(STR) == EXPECTED
        where:
        STR         | EXPECTED
        'abc'       | 'abc'
        'a_b_c'     | 'aBC'
        'foo__bar'  | 'fooBar'
    }
    

    def 'should validate URL' () {
        given:
        def observer = new TowerClient()

        expect:
        observer.checkUrl('http://localhost') == 'http://localhost'
        observer.checkUrl('http://google.com') == 'http://google.com'
        observer.checkUrl('https://google.com') == 'https://google.com'
        observer.checkUrl('http://google.com:8080') == 'http://google.com:8080'
        observer.checkUrl('http://google.com:8080/') == 'http://google.com:8080'
        observer.checkUrl('http://google.com:8080/foo/bar') == 'http://google.com:8080/foo/bar'
        observer.checkUrl('http://google.com:8080/foo/bar/') == 'http://google.com:8080/foo/bar'
        observer.checkUrl('http://google.com:8080/foo/bar///') == 'http://google.com:8080/foo/bar'

        when:
        observer.checkUrl('ftp://localhost')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Only http and https are supported -- The given URL was: ftp://localhost'
    }

    def 'should get watch url' () {
        given:
        def observer = new TowerClient()
        expect:
        observer.getHostUrl(STR) == EXPECTED
        where:
        STR                             | EXPECTED
        'http://foo.com'                | 'http://foo.com'
        'http://foo.com:800/'           | 'http://foo.com:800'
        'https://foo.com:800/'          | 'https://foo.com:800'
        'http://foo.com:8000/this/that' | 'http://foo.com:8000'
    }

    def 'should get access token' () {
        given:
        def session = Mock(Session)

        when:
        def config = new TowerConfig([accessToken: 'abc'], [TOWER_ACCESS_TOKEN: 'xyz'])
        def observer = new TowerClient(session, config)
        then:
        // the token in the config overrides the one in the env
        observer.getAccessToken() == 'abc'

        when:
        config = new TowerConfig([accessToken: 'abc'], [TOWER_ACCESS_TOKEN: 'xyz', TOWER_WORKFLOW_ID: '111222333'])
        observer = new TowerClient(session, config)
        then:
        // the token from the env is taken because is a tower launch aka TOWER_WORKFLOW_ID is set
        observer.getAccessToken() == 'xyz'

        when:
        config = new TowerConfig([:], [TOWER_ACCESS_TOKEN: 'xyz'])
        observer = new TowerClient(session, config)
        then:
        observer.getAccessToken() == 'xyz'

        when:
        config = new TowerConfig([:], [:])
        observer = new TowerClient(session, config)
        observer.getAccessToken()
        then:
        thrown(AbortOperationException)
    }

    def 'should post task records' () {
        given:
        def URL = 'http://foo.com'
        def PROGRESS = Mock(WorkflowProgress) { getRunning()>>1; getSucceeded()>>2; getFailed()>>3 }
        def client = Mock(HxClient)
        def observer = Spy(TowerClient)
        observer.@httpClient = client
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

        when:
        observer.sendHttpMessage(URL, req)
        then:
        1 * client.sendAsString(_) >> Mock(HttpResponse)

    }


    static now_millis = System.currentTimeMillis()
    static now_instant = OffsetDateTime.ofInstant(Instant.ofEpochMilli(now_millis), ZoneId.systemDefault())

    def 'should fix field types' () {

        expect:
        TowerClient.fixTaskField(FIELD,VALUE) == EXPECTED

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
        def http = Mock(HxClient)
        TowerClient client = Spy(new TowerClient([httpClient: http, env: ENV]))
        and:
        client.getOperationId() >> 'op-112233'
        client.getLogFile() >> 'log.file'
        client.getOutFile() >> 'out.file'

        and:
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

        when:
        def req1 = client.makeCreateReq(session)
        then:
        req1.sessionId == sessionId.toString()
        req1.runName == 'foo'
        req1.projectName == 'the-project-name'
        req1.repository == 'git://repo.com/foo'
        req1.workflowId == WORKFLOW_ID
        and:
        aroundNow(req1.instant)

        when:
        def req = client.makeBeginReq(session)
        then:
        client.getWorkflowId() >> '12345'
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
        def tower = new TowerClient()

        expect:
        tower.mapToString(null)  == null
        tower.mapToString('ciao') == 'ciao'
        tower.mapToString([:]) == null
        tower.mapToString([p:'foo', q:'bar']) == null
    }


    def 'should load schema col len' () {
        given:
        def tower = new TowerClient()

        when:
        def schema = tower.loadSchema()
        then:
        schema.get('workflow.start')  == null
        schema.get('workflow.profile') == 100
        schema.get('workflow.projectDir') == 255
    }

    def 'should create init request' () {
        given:
        def uuid = UUID.randomUUID()
        def client = new TowerClient(env: [TOWER_WORKFLOW_ID: 'x123'])
        def meta = Mock(WorkflowMetadata) {
            getProjectName() >> 'the-project-name'
            getRepository() >> 'git://repo.com/foo'
        }
        def session = Mock(Session) {
            getUniqueId() >> uuid
            getRunName() >> 'foo_bar'
            getWorkflowMetadata() >> meta
        }

        when:
        def req = client.makeCreateReq(session)
        then:
        req.sessionId == uuid.toString()
        req.runName == 'foo_bar'
        req.projectName == 'the-project-name'
        req.repository == 'git://repo.com/foo'
        req.workflowId == 'x123'
        and:
        aroundNow(req.instant)

        and:
        client.towerLaunch
    }

    def 'should post create request' () {
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
        def config = new TowerConfig([:], [:])

        def client = Spy(new TowerClient(session, config))

        when:
        client.onFlowCreate(session)
        then:
        1 * client.getAccessToken() >> 'secret'
        1 * client.makeCreateReq(session) >> [runName: 'foo']
        1 * client.sendHttpMessage('https://api.cloud.seqera.io/trace/create', [runName: 'foo'], 'POST') >> new TowerClient.Response(200, '{"workflowId":"xyz123"}')
        and:
        client.runName == 'foo_bar'
        client.runId == uuid.toString()
        and:
        client.workflowId == 'xyz123'
        !client.towerLaunch

    }

    def 'should get trace endpoint' () {
        given:
        def config = new TowerConfig([:], [:])
        def tower = new TowerClient(Mock(Session), config)
        tower.workflowId = '12345'

        expect:
        tower.getUrlTraceCreate() == 'https://api.cloud.seqera.io/trace/create'
        tower.getUrlTraceBegin() == 'https://api.cloud.seqera.io/trace/12345/begin'
        tower.getUrlTraceProgress() == 'https://api.cloud.seqera.io/trace/12345/progress'
        tower.getUrlTraceHeartbeat() == 'https://api.cloud.seqera.io/trace/12345/heartbeat'
        tower.getUrlTraceComplete() == 'https://api.cloud.seqera.io/trace/12345/complete'
    }

    def 'should get trace endpoint with workspace' () {
        given:
        def config = new TowerConfig([workspaceId: '300'], [:])
        def tower = new TowerClient(Mock(Session), config)
        tower.workflowId = '12345'

        expect:
        tower.getUrlTraceCreate() == 'https://api.cloud.seqera.io/trace/create?workspaceId=300'
        tower.getUrlTraceBegin() == 'https://api.cloud.seqera.io/trace/12345/begin?workspaceId=300'
        tower.getUrlTraceProgress() == 'https://api.cloud.seqera.io/trace/12345/progress?workspaceId=300'
        tower.getUrlTraceHeartbeat() == 'https://api.cloud.seqera.io/trace/12345/heartbeat?workspaceId=300'
        tower.getUrlTraceComplete() == 'https://api.cloud.seqera.io/trace/12345/complete?workspaceId=300'
    }

    def 'should set the auth token' () {
        given:
        def http = Mock(HxClient.Builder)
        def session = Mock(Session)
        def config = new TowerConfig([:], [:])
        def client = Spy(new TowerClient(session, config))
        and:
        def SIMPLE = '4ffbf1009ebabea77db3d72efefa836dfbb71271'
        def BEARER = 'eyJ0aWQiOiA1fS5jZmM1YjVhOThjZjM2MTk1NjBjZWU1YmMwODUxYzA1ZjkzMDdmN2Iz'

        when:
        client.setupClientAuth(http, SIMPLE)
        then:
        1 * http.basicAuth('@token:' + SIMPLE) >> http

        when:
        client.setupClientAuth(http, SIMPLE)
        then:
        1 * http.basicAuth('@token:' + SIMPLE) >> http

        when:
        client.setupClientAuth(http, BEARER)
        then:
        1 * http.bearerToken(BEARER) >> http
        1 * http.refreshToken(_) >> http
        1 * http.refreshTokenUrl(_) >> http
    }

    def 'should fetch workflow meta' () {
        given:
        def client = Spy(new TowerClient(env: ENV))

        expect:
        client.getOperationId() == OP_ID
        client.getLogFile() == LOG_FILE
        client.getOutFile() == OUT_FILE

        where:
        OP_ID                                           | OUT_FILE      | LOG_FILE    | ENV
        null                                            | null          | null        | [:]
        "local-platform::${ProcessHelper.selfPid()}"    | null          | null        | [TOWER_ALLOW_NEXTFLOW_LOGS:'true']
        'aws-batch::1234z'                              | 'xyz.out'     | 'hola.log'  | [TOWER_ALLOW_NEXTFLOW_LOGS:'true', AWS_BATCH_JOB_ID: '1234z', NXF_OUT_FILE: 'xyz.out', NXF_LOG_FILE: 'hola.log']
    }

    def 'should deduplicate containers' () {
        given:
        def client = Spy(new TowerClient())
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
        client.getNewContainers([trace1]) == [c1]
        and:
        client.getNewContainers([trace1]) == []
        and:
        client.getNewContainers([trace1, trace2, trace3]) == [c2]
    }

    def 'should handle HTTP request with content'() {
        given: 'a TowerClient'
        def tower = new TowerClient()
        def content = '{"test": "data"}'
        def request = tower.makeRequest('http://example.com/test', content, 'POST')

        expect: 'the request should be created with the content'
        request != null
        request.method() == 'POST'
        request.uri().toString() == 'http://example.com/test'
    }

    def 'should include numSpotInterruptions in task map'() {
        given:
        def client = Spy(new TowerClient())
        client.getWorkflowProgress(true) >> new WorkflowProgress()

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
        trace.setnumSpotInterruptions(3)

        when:
        def req = client.makeTasksReq([trace])

        then:
        req.tasks.size() == 1
        req.tasks[0].numSpotInterruptions == 3
    }

}
