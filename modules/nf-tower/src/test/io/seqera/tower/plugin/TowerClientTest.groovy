/*
 * Copyright (c) 2019, Seqera Labs.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */

package io.seqera.tower.plugin

import java.nio.file.Files
import java.nio.file.Paths
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneId

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.exception.AbortOperationException
import nextflow.executor.LocalExecutor
import nextflow.executor.LocalTaskHandler
import nextflow.processor.TaskConfig
import nextflow.processor.TaskId
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.WorkflowMetadata
import nextflow.util.SimpleHttpClient
import org.junit.Ignore
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerClientTest extends Specification {

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
        tower.runName = session.runName
        tower.workflowId = '12ef'
        tower.terminated = true

        when:
        def map = tower.makeWorkflowReq(session)
        then:
        1 * session.getWorkflowMetadata() >> meta
        1 * session.getParams() >> params
        1 * meta.toMap() >> [foo:1, bar:2, container: [p1: 'c1', p2: 'c2']]
        1 * tower.getMetricsList() >> [[process:'foo', cpu: [min: 1, max:5], time: [min: 6, max: 9]]]
        then:
        map.workflow.foo == 1
        map.workflow.bar == 2
        map.workflow.id == '12ef'
        map.workflow.params == [x: 'hello']
        map.workflow.container == 'p1:c1,p2:c2'
        map.metrics == [[process:'foo', cpu: [min: 1, max:5], time: [min: 6, max: 9]]]

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
        e.message == 'Only http or https are supported protocols -- The given URL was: ftp://localhost'
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
        def ENV = [TOWER_ACCESS_TOKEN: 'xyz']
        def session = Mock(Session)

        when:
        def observer = new TowerClient(session: session)
        def result = observer.getAccessToken()
        then:
        session.getConfig() >> [tower:[accessToken: 'abc'] ]
        and:
        result == 'abc'

        when:
        observer = new TowerClient(session: session, env: ENV)
        result = observer.getAccessToken()
        then:
        session.getConfig() >> [:]
        and:
        result == 'xyz'

        when:
        observer = new TowerClient(session: session, env:[:])
        observer.getAccessToken()
        then:
        session.getConfig() >> [:]
        then:
        thrown(AbortOperationException)

    }

    def 'should post task records' () {
        given:
        def URL = 'http://foo.com'
        def client = Mock(SimpleHttpClient)
        def session = Mock(Session)
        def proc = Mock(TaskProcessor) { getName() >> 'foo'; getProcessEnvironment() >> [:] }
        def observer = Spy(TowerClient)
        observer.session = session
        observer.httpClient = client

        def now = System.currentTimeMillis()

        def task = new TaskRun(
                id: TaskId.of(10),
                workDir: Paths.get('/work/dir'),
                config: new TaskConfig(),
                processor: proc
        )
        def handler = new LocalTaskHandler(task, Mock(LocalExecutor))
        handler.submitTimeMillis = now-2000
        handler.startTimeMillis = now-1000
        handler.completeTimeMillis = now
        def trace = handler.getTraceRecord()

        when:
        def req = observer.makeTasksReq([trace])
        then:
        req.tasks[0].taskId == 10
        req.tasks[0].process == 'foo'
        req.tasks[0].workdir == "/work/dir"
        req.tasks[0].cpus == 1
        req.tasks[0].submit == OffsetDateTime.ofInstant(Instant.ofEpochMilli(handler.submitTimeMillis), ZoneId.systemDefault())
        req.tasks[0].start == OffsetDateTime.ofInstant(Instant.ofEpochMilli(handler.startTimeMillis), ZoneId.systemDefault())

        when:
        observer.sendHttpMessage(URL, req)
        then:
        1 * client.sendHttpMessage(URL, _, 'POST') >> null

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


    @Ignore
    def 'should create workflow json' () {

        given:
        def dir = Files.createTempDirectory('test')
        def client = Mock(SimpleHttpClient)
        TowerClient tower = Spy(TowerClient, constructorArgs: [[httpClient: client] ])

        def session = Mock(Session)
        session.config >> [:]
        session.containerConfig >> new ContainerConfig()
        session.getParams() >> new ScriptBinding.ParamsMap([foo:'Hello', bar:'World'])

        def meta = new WorkflowMetadata(session, Mock(ScriptFile))
        session.getWorkflowMetadata() >> meta
        when:
        def req = tower.makeWorkflowReq(session)
        println req.workflow
        then:
        req.workflow.params == [foo:'Hello', bar:'World']


        cleanup:
        dir?.deleteDir()

    }

    def 'should convert map' () {
        given:
        def tower = new TowerClient()

        expect:
        tower.mapToString(null)  == null
        tower.mapToString('ciao') == 'ciao'
        tower.mapToString([:]) == ''
        tower.mapToString([p:'foo', q:'bar']) == 'p:foo,q:bar'
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
        def tower = new TowerClient()
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
        def req = tower.makeInitRequest(session)
        then:
        req.sessionId == uuid.toString()
        req.runName == 'foo_bar'
        req.projectName == 'the-project-name'
        req.repository == 'git://repo.com/foo'
    }

    def 'should post init request' () {
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

        TowerClient tower = Spy(TowerClient, constructorArgs: ['https://tower.nf'])

        when:
        tower.onFlowInit(session)
        then:
        1 * tower.getAccessToken() >> 'secret'
        1 * tower.makeInitRequest(session) >> [runName: 'foo']
        1 * tower.sendHttpMessage('https://tower.nf/trace/init', [runName: 'foo']) >> new TowerClient.Response(200, '{"workflowId":"xyz123"}')
        and:
        tower.runName == 'foo_bar'
        tower.runId == uuid.toString()
        and:
        tower.workflowId == 'xyz123'

    }
}
