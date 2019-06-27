/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.processor

import java.util.concurrent.ArrayBlockingQueue
import org.codehaus.groovy.runtime.powerassert.PowerAssertionError

import nextflow.Session
import nextflow.util.Duration
import nextflow.util.RateUnit
import nextflow.wr.client.WrRestApi
import nextflow.wr.executor.WrTaskHandler
import spock.lang.Specification
/**
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TaskPollingMonitorTest Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrMonitorTest extends Specification {

    def 'should create a wr monitor'() {
        setup:
        def session = new Session( executor: [pollInterval: '1h', queueSize: 11, dumpInterval: '3h'] )
        def client = Mock(WrRestApi)

        when:
        def monitor = WrMonitor.create(session, client)

        then:
        monitor.client == client
        monitor.name == 'wr'
        monitor.capacity == 11
        monitor.pollIntervalMillis == Duration.of('1h').toMillis()
        monitor.dumpInterval ==  Duration.of('3h')

        when:
        monitor = new WrMonitor(name:'local', session: session, pollInterval: '1s', capacity: 100)

        then:
        thrown(PowerAssertionError)
    }

    def 'should submit'() {
        given:
        def session = Spy(new Session( executor: [pollInterval: '1h', queueSize: 11, dumpInterval: '3h'] ))
        def client = Mock(WrRestApi)
        def monitor = Spy(WrMonitor.create(session, client))
        def handler = Mock(WrTaskHandler)
        String id = "foo"

        when:
        monitor.submit(handler, id)

        then:
        1 * handler.submitted(id)
        // 1 * monitor.runningQueue.add(handler) *** don't know how to test this happens
        1 * session.notifyTaskSubmit(handler) >> null
    }

    def 'should return true for canSubmit' () {
        setup:
        def session = new Session( executor: [pollInterval: '1h', queueSize: 11, dumpInterval: '3h'] )
        def client = Mock(WrRestApi)
        def monitor = WrMonitor.create(session, client)
        def handler = Mock(WrTaskHandler)

        when:
        def can = monitor.canSubmit(handler)

        then:
        can == true
    }

    def 'should submit all pending tasks' () {
        setup:
        def session = Spy(new Session( executor: [pollInterval: '1h', queueSize: 11, dumpInterval: '3h'] ))
        List<List> toSubmit = [["exe one"], ["exe two"]]
        List<Map> jobs = []
        jobs << ["Cmd":toSubmit[0][0], "Key":"key1"]
        jobs << ["Cmd":toSubmit[1][0], "Key":"key2"]
        WrRestApi client = Mock {
            1 * add(toSubmit) >> jobs
        }
        def monitor = Spy(WrMonitor.create(session, client))
        WrTaskHandler handler1 = Mock {
            1 * submitArgs() >> toSubmit[0]
        }
        WrTaskHandler handler2 = Mock {
            1 * submitArgs() >> toSubmit[1]
        }

        when:
        monitor.start()
        monitor.schedule(handler1)
        monitor.schedule(handler2)
        def submitted = monitor.submitPendingTasks()

        then:
        submitted == 2
        1 * monitor.pollLoop() >> null
        1 * monitor.submitLoop() >> null
        1 * monitor.submit(handler1, "key1")
        1 * monitor.submit(handler2, "key2")
        1 * session.notifyTaskSubmit(handler1) >> null
        1 * session.notifyTaskSubmit(handler2) >> null
    }

}
