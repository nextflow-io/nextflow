/*
 * Copyright (c) 2018, University of Tübingen, Quantitative Biology Center (QBiC).
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.trace

import spock.lang.Specification

import groovy.json.JsonSlurper
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.util.SimpleHttpClient

class WebLogObserverTest extends Specification{

    def 'do not send messages on wrong formatted url' () {

        when:
        new WebLogObserver("localhost")

        then:
        thrown(IllegalArgumentException)
    }

    def 'send message on different workflow events' () {

        given:
        def httpPostObserver0 = Spy(WebLogObserver, constructorArgs: ["http://localhost"])
        def sessionStub = Mock(Session){
            getRunName() >> "testRun"
            getUniqueId() >> UUID.randomUUID()
        }
        def traceStub = Mock(TraceRecord)
        def handlerStub = Mock(TaskHandler)

        when:
        httpPostObserver0.onFlowStart(sessionStub)
        httpPostObserver0.onFlowComplete()
        httpPostObserver0.onProcessSubmit(handlerStub, traceStub)
        httpPostObserver0.onProcessStart(handlerStub, traceStub)
        httpPostObserver0.onProcessComplete(handlerStub, traceStub)
        httpPostObserver0.onFlowError(handlerStub, traceStub)

        then:
        noExceptionThrown()
        2 * httpPostObserver0.asyncHttpMessage(_)
        4 * httpPostObserver0.asyncHttpMessage(!null, !null)

    }

    def 'should create a json message' () {

        given:
        def observer = Spy(WebLogObserver)
        def CLIENT = Mock(SimpleHttpClient)
        observer.httpClient = CLIENT
        observer.runName = 'foo'
        observer.runId = 'xyz'
        def TRACE = new TraceRecord([hash: '4a4a4a', process: 'bar'])

        when:
        observer.sendHttpMessage('started',  TRACE)
        then:
        1 * observer.logHttpResponse() >> null
        1 * CLIENT.sendHttpMessage( _ as String ) >> { it ->
            def message = (Map)new JsonSlurper().parseText((String)it[0])
            assert message.runName == 'foo'
            assert message.runId == 'xyz'
            assert message.event == 'started'
            assert message.runStatus == 'started'
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
        e.message == 'Only http or https are supported protocols -- The given URL was: ftp://localhost'
    }
}
