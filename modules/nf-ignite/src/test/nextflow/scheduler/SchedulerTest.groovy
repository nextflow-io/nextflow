/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.scheduler

import static nextflow.scheduler.Protocol.PENDING_TASKS_CACHE
import static nextflow.scheduler.Protocol.TOPIC_AGENT_EVENTS

import nextflow.executor.IgBaseTask
import nextflow.processor.TaskId
import nextflow.processor.TaskPollingMonitor
import nextflow.util.Duration
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteCompute
import org.apache.ignite.IgniteEvents
import org.apache.ignite.IgniteMessaging
import org.apache.ignite.events.EventType
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SchedulerTest extends Specification {


    def 'should send a task for execution' () {

        given:
        def now = System.currentTimeMillis()
        def _10_mins_ago = now - Duration.of('10 min').toMillis()
        def _57_mins_ago = now - Duration.of('57 min').toMillis()
        def _90_mins_ago = now - Duration.of('90 min').toMillis()
        def n1 = UUID.fromString('11111111-9b9f-4b45-bebc-aae51fad9d18')
        def n2 = UUID.fromString('22222222-2b3b-4a98-b771-f20d1ed50c83')
        def n3 = UUID.fromString('33333333-663c-4760-8300-4d0ff90ef7cd')
        def n4 = UUID.fromString('44444444-663c-4760-8300-4d0ff90ef7cd')
        def n5 = UUID.fromString('55555555-663c-4760-8300-4d0ff90ef7cd')

        def topology = [
                new Protocol.NodeData(nodeId: n1, instanceId: 'i1', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10),
                new Protocol.NodeData(nodeId: n2, instanceId: 'i2', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10_000),
                new Protocol.NodeData(nodeId: n3, instanceId: 'i3', bootTimeMillis: _10_mins_ago, idleTimestamp: now-10_000),
                new Protocol.NodeData(nodeId: n4, instanceId: 'i4', bootTimeMillis: _57_mins_ago, idleTimestamp: now-10_000),
                new Protocol.NodeData(nodeId: n5, instanceId: 'i5', bootTimeMillis: _90_mins_ago, idleTimestamp: now-10_000),
        ]


        def scheduler = Spy(Scheduler)

        def messaging = Mock(IgniteMessaging)
        def events = Mock(IgniteEvents)
        def monitor = Mock(TaskPollingMonitor)
        def ignite = Mock(Ignite)
        def compute = Mock(IgniteCompute)
        def pendingTasks = Mock(IgniteCache)
        ignite.message() >> messaging
        ignite.events() >> events
        ignite.compute() >> compute

        def task = Mock(IgBaseTask)
        task.getTaskId() >> TaskId.of(11)

        when:
        scheduler.init( ignite, monitor )
        then:
        1 * ignite.cache(PENDING_TASKS_CACHE) >> pendingTasks
        1 * compute.broadcast( _ ) >> topology
        1 * messaging.localListen(Protocol.TOPIC_SCHEDULER_EVENTS, _)
        1 * events.localListen(_, EventType.EVT_NODE_LEFT)
        1 * events.localListen(_, EventType.EVT_NODE_FAILED)
        scheduler.isRunning()

        when:
        scheduler.schedule0(task)
        then:
        1 * pendingTasks.put(task.getTaskId(), task)
        1 * messaging.send(TOPIC_AGENT_EVENTS, Protocol.TaskAvail.INSTANCE)
    }

}
