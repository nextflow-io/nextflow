/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2018, University of Tübingen, Quantitative Biology Center (QBiC)
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.trace

import spock.lang.Specification

import groovy.json.JsonGenerator
import groovy.json.JsonSlurper
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.script.ScriptBinding
import nextflow.script.WorkflowMetadata
import nextflow.util.SimpleHttpClient

class WebLogObserverTest extends Specification {

    def 'do not send messages on wrong formatted url'() {

        when:
        new WebLogObserver("localhost")

        then:
        thrown(IllegalArgumentException)
    }

    def 'send message on different workflow events' () {

        given:
        WebLogObserver httpPostObserver0 = Spy(WebLogObserver, constructorArgs: ["http://localhost"])
        WorkflowMetadata workflowMeta = Mock(WorkflowMetadata)

        def bindingStub = Mock(ScriptBinding){
            getProperty('params') >> new ScriptBinding.ParamsMap()
        }
        def sessionStub = Mock(Session){
            getRunName() >> "testRun"
            getUniqueId() >> UUID.randomUUID()
            getBinding() >> bindingStub
            getWorkflowMetadata() >> workflowMeta
        }
        def traceStub = Mock(TraceRecord)
        def handlerStub = Mock(TaskHandler)

        when:
        def payload = WebLogObserver.createFlowPayloadFromSession(sessionStub)
        httpPostObserver0.onFlowCreate(sessionStub)
        httpPostObserver0.onProcessSubmit(handlerStub, traceStub)
        httpPostObserver0.onProcessStart(handlerStub, traceStub)
        httpPostObserver0.onProcessComplete(handlerStub, traceStub)
        httpPostObserver0.onFlowError(handlerStub, traceStub)
        httpPostObserver0.onFlowComplete()

        then:
        assert payload.workflow instanceof WorkflowMetadata
        6 * httpPostObserver0.asyncHttpMessage(!null, !null) >> null

    }

    def 'should create a json message' () {

        given:
        def observer = Spy(WebLogObserver)
        def CLIENT = Mock(SimpleHttpClient)
        observer.@endpoint = 'http://foo.com'
        observer.@httpClient = CLIENT
        observer.@runName = 'foo'
        observer.@runId = 'xyz'
        observer.@generator = new JsonGenerator.Options().build()
        def TRACE = new TraceRecord([hash: '4a4a4a', process: 'bar'])

        when:
        observer.sendHttpMessage('started', TRACE)

        then:
        1 * observer.logHttpResponse() >> null
        1 * CLIENT.sendHttpMessage( 'http://foo.com', _ as String ) >> { it ->
            def message = (Map)new JsonSlurper().parseText((String)it[1])
            assert message.runName == 'foo'
            assert message.runId == 'xyz'
            assert message.event == 'started'
            assert message.trace.hash == '4a4a4a'
            assert message.trace.process == 'bar'
            return null
        }

    }

    def 'should validate URL' () {
        given:
        def observer = new WebLogObserver()
        
        expect:
        observer.checkUrl('http://localhost') == 'http://localhost'
        observer.checkUrl('http://google.com') == 'http://google.com'
        observer.checkUrl('https://google.com') == 'https://google.com'
        observer.checkUrl('http://google.com:8080') == 'http://google.com:8080'
        observer.checkUrl('http://google.com:8080/foo/bar') == 'http://google.com:8080/foo/bar'

        when:
        observer.checkUrl('ftp://localhost')
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Only http and https are supported -- The given URL was: ftp://localhost'
    }
}
