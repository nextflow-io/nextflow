/*
 * Copyright (c) 2012, the authors.
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

package nextflow.executor

import java.util.concurrent.Callable

import com.hazelcast.client.HazelcastClient
import com.hazelcast.client.config.ClientConfig
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IExecutorService
import com.hazelcast.core.IQueue
import com.hazelcast.core.ITopic
import com.hazelcast.core.Member
import com.hazelcast.core.MembershipEvent
import com.hazelcast.core.MembershipListener
import com.hazelcast.core.Message
import com.hazelcast.core.MessageListener
import com.hazelcast.core.MultiExecutionCallback
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskPollingMonitor

/**
 * Creates a connector for the Hazelcast cluster
 *
 * @author Paolo Di Tommaso
 */
@Slf4j
class HzConnector implements HzConst, MembershipListener, MessageListener<HzCmdResult> {

    private HazelcastInstance hazelcast

    private ITopic<HzCmdResult> resultsTopic

    private IQueue<HzCmdCall> executorsQueue

    @PackageScope
    private final Map<Member,Integer> slotsFor = [:]


    private TaskPollingMonitor monitor

    /** The nextflow session object */
    private Session session

    private Map<UUID,HzRemoteSession> allSessions


    /**
     * Factory method. Use this method to create an instance of the Hazelcast connect
     */
    @Memoized
    static HzConnector create( TaskPollingMonitor monitor ) {
        log.debug "Creating Hazelcast connector object for monitor: $monitor"
        new  HzConnector( monitor, monitor.session)
    }

    /**
     * Class constructor
     *
     * @param monitor
     * @param session
     */
    protected HzConnector( TaskPollingMonitor monitor, Session session ) {
        this.monitor = monitor
        this.session = session
        initialize()

    }


    TaskPollingMonitor getTaskMonitor() {
        monitor
    }


    @PackageScope
    Queue<HzCmdCall> getExecutorsQueue() {
        return executorsQueue
    }


    /**
     * Initialize the Hazelcast client connection
     *
     * @return A {@code HazelcastInstance} instance
     */
    protected void initialize() {
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')

        def cfg = new ClientConfig()
        hazelcast = HazelcastClient.newHazelcastClient(cfg)

        executorsQueue = hazelcast.getQueue(TASK_SUBMITS_NAME)
        resultsTopic = hazelcast.getTopic(TASK_RESULTS_NAME)
        resultsTopic.addMessageListener(this)
        // publish the session on the cluster
        allSessions = hazelcast.getMap(SESSION_MAP)
        allSessions.put(session.uniqueId, new HzRemoteSession(session)  )

        // add the members listener
        final cluster = hazelcast.getCluster()
        cluster.addMembershipListener( this )
        Set members  = cluster.getMembers();
        log.debug "Cluster members (${members.size()}): ${ members.size() <= 10 ? members.join(', ') : members.take(10).join(', ') + " ..." }"

        // fetch the current number of slots for each member
        remote().submitToAllMembers(new AvailCores(), new AvailCoresCallback(this))

        // register for the shutdown
        session.onShutdown {
            if( !hazelcast ) return
            log.info "Shutting down hazelcast connector"
            allSessions.remove(session.getUniqueId())
            hazelcast.shutdown()
        }
    }


    /**
     * This method is invoked for message on topic {@code #TASK_RESULTS_NAME}
     *
     * @param message The {@code HzCmdResult} represent the result of a task executed in the Hazelcast cluster
     */
    @Override
    void onMessage(Message<HzCmdResult> message) {
        log.trace "Received command message: $message"

        // -- make sure the sender match with the current session id
        def result = message.getMessageObject()
        if( result.sender != session.uniqueId ) {
            log.debug "Discarding command result: $result"
            return
        }

        // -- find out the task handler for the task id of the completed task
        log.trace "Received command result: $result"
        def handler = (HzTaskHandler)taskMonitor.getTaskHandlerBy(result.taskId)

        if( !handler ) {
            log.warn "Lost task for command result: $result"
            return
        }

        // -- set the result value and *signal* the monitor in order to trigger a status check cycle
        log.trace "Setting result: $result to handler: $handler"
        handler.result = result
        taskMonitor.signalComplete()

    }

    /**
     * Invoked by Hazelcast when a new cluster member is added
     *
     * @param membershipEvent
     */
    @Override
    void memberAdded(MembershipEvent e) {
        log.debug "Cluster member added: ${e.member}"

        final result = remote().submitToMember(new AvailCores(), e.member)
        final slots = result.get()

        log.debug "Cluster member: ${e.member} provides ${slots} slots"
        slotsFor[e.member] = slots

        // -- increase the monitor capacity by the number of slots provided by the member added
        if( slots )
            taskMonitor.capacityInc(slots)

    }

    /**
     * Invoked by Hazelcast when a cluster member is removed
     *
     * @param membershipEvent
     */
    @Override
    void memberRemoved(MembershipEvent e) {
        log.debug "Cluster member removed: ${e.member}"

        final slots = slotsFor[e.member]
        if( !slots ) {
            log.debug "Unknown slots for member: ${e.member} -- won't reduce capacity"
            return
        }

        // -- reduce the monitor capacity by the number of slots provided by the member removed
        slotsFor.remove(e.member)
        taskMonitor.capacityDec(slots)

    }


    @PackageScope
    IExecutorService remote() {
        hazelcast.getExecutorService(EXEC_SERVICE)
    }



    /**
     * A callable to get the number of cores/slots from a remote Hazelcast node
     */
    static class AvailCores implements Callable<Integer>, Serializable {

        private static final long serialVersionUID = - 7880952114196721598L ;

        @Override
        Integer call() throws Exception {
            Runtime.getRuntime().availableProcessors()
        }
    }

    /**
     * Notifies when the {@code AvailCores} remote call has completed
     */
    @Slf4j
    static class AvailCoresCallback implements MultiExecutionCallback {

        private TaskPollingMonitor monitor

        private HzConnector connector

        AvailCoresCallback( HzConnector connector ) {
            assert connector
            this.connector = connector
            this.monitor = connector.getTaskMonitor()
        }

        @Override
        void onResponse(Member member, Object value) {
            if( member && value instanceof Integer ) {
                log.debug "Increasing hazelcast monitor capacity of $value slots provided by $member"
                connector.slotsFor[member] = value
                monitor.capacityInc(value as int)
            }

            else
                log.debug "Invalid AvailCoresCallback data > member: $member; value: $value"
        }

        @Override
        void onComplete(Map<Member, Object> values) {
            // ignore
        }
    }

}