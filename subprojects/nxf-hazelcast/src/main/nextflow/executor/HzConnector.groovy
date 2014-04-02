/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import com.hazelcast.client.HazelcastClient
import com.hazelcast.client.config.ClientConfig
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IMap
import com.hazelcast.core.IQueue
import com.hazelcast.core.ItemEvent
import com.hazelcast.core.ItemEventType
import com.hazelcast.core.ItemListener
import com.hazelcast.core.Member
import com.hazelcast.core.MemberAttributeEvent
import com.hazelcast.core.MembershipEvent
import com.hazelcast.core.MembershipListener
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskPollingMonitor
import org.apache.commons.lang.StringUtils
/**
 * Creates a connector for the Hazelcast cluster
 *
 * @author Paolo Di Tommaso
 */
@Slf4j
class HzConnector implements HzConst, MembershipListener, ItemListener {

    private HazelcastInstance hazelcast

    private IQueue<HzCmdCall> submitQueue

    private IQueue<HzCmdNotify> eventsQueue

    private IMap<String,HzNodeInfo> allMembers

    private TaskPollingMonitor monitor

    /** The nextflow session object */
    private Session session

    private Map<UUID,HzRemoteSession> allSessions

    private capacityFactory = 1.5f


    /**
     * Factory method. Use this method to create an instance of the Hazelcast connect
     */
    @Memoized
    static HzConnector create( TaskPollingMonitor monitor ) {
        log.debug "Creating Hazelcast connector object for monitor: $monitor"
        new  HzConnector( monitor, monitor.session )
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

    /** Only for test -- DO NOT USE */
    protected HzConnector() {}


    @PackageScope
    Queue<HzCmdCall> getExecutorsQueue() {
        return submitQueue
    }

    protected getHzConfigProperty( String key, defValue = null )  {
        session.getExecConfigProp('hazelcast', key, defValue)
    }

    protected String getHzGroup() {
        getHzConfigProperty('group', DEFAULT_GROUP_NAME)
    }

    protected List<String> getHzJoinAddress() {
        def result = getHzConfigProperty('join')
        return result ? StringUtils.split(result.toString(), ", \n").collect { it.trim()  } : []
    }

    /**
     * Initialize the Hazelcast client connection
     *
     * @return A {@code HazelcastInstance} instance
     */
    protected void initialize() {
        // -- configure the logging
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')

        // -- get info from the config file available from the session
        def cfg = new ClientConfig()
        def group = getHzGroup()
        log.debug "Hz client config > group: $group"
        cfg.getGroupConfig().setName(group)
        getHzJoinAddress().each {
            log.debug "Hz client config > Adding join address: $it"
            cfg.addAddress(it)
        }

        /*
         * create the main Hazelcast instance
         */
        HzSerializerConfig.registerAll(cfg.getSerializationConfig())
        hazelcast = HazelcastClient.newHazelcastClient(cfg)

        // -- queue where tasks are submitted
        submitQueue = hazelcast.getQueue(TASK_SUBMITS_QUEUE)

        eventsQueue = hazelcast.getQueue(TASK_EVENTS_QUEUE)
        eventsQueue.addItemListener(this,true)

        // publish the session on the cluster
        allSessions = hazelcast.getMap(SESSIONS_MAP)
        allSessions.put(session.uniqueId, new HzRemoteSession(session)  )

        allMembers = hazelcast.getMap(MEMBERS_MAP)
        def capacity = 0
        allMembers.values().each { HzNodeInfo it -> capacity += it.slots }
        log.debug "Cluster initial capacity > $capacity"
        capacityInc(capacity)

        // add the members listener
        final cluster = hazelcast.getCluster()
        cluster.addMembershipListener( this )
        Set members  = cluster.getMembers();
        log.debug "Cluster members (${members.size()}): ${ members.size() <= 10 ? members.join(', ') : members.take(10).join(', ') + " ..." }"

        // register for the shutdown
        session.onShutdown {
            if( !hazelcast ) return
            log.info "Shutting down hazelcast connector"
            allSessions.remove(session.getUniqueId())
            hazelcast.shutdown()
        }
    }

    private void capacityInc( int value ) {
        log.debug "Incrementing capacity by $value slots"
        final capacity = (value * capacityFactory as int) +1
        monitor.capacityInc(capacity)
    }

    private void capacityDec( int value ) {
        log.debug "Decremeting capacity by $value slots"
        final capacity = (value * capacityFactory as int) +1
        monitor.capacityDec(capacity)
    }

    /*
     * Get invoked when a new task is queued by a remote daemon
     */
    @Override
    void itemAdded(ItemEvent event) {

        if( event?.eventType != ItemEventType.ADDED ) {
            log.debug "Not interested on this event: $event"
        }

        HzCmdNotify message = (HzCmdNotify)event.item
        log.trace "Received notification: $message"
        if( !message ) {
            return
        }

        if( message.sessionId != session.uniqueId ) {
            log.trace "Discarding command notification since belongs to other client: $message"
            return
        }

        // just consume this entry
        eventsQueue.remove(message)

        /*
         * find out the task handler for the task id of the completed task
         */
        def handler = (HzTaskHandler)monitor.getTaskHandlerBy(message.taskId)
        if( !handler ) {
            log.warn "Lost task for command start: ${message.taskId}"
            return
        }

        /*
         * when it is a task completion notification,
         * update the handler with that object
         */
        if( message.isComplete() ) {
            log.trace "Task complete > ${message}"
            monitor.schedule {
                handler.result = message
            }
            return
        }

        /*
         * .. otherwise it is an task task notification
         */

        /*
         * make sure the cluster member is still available
         */
        if( isAvailable(message.memberId) ) {
            // -- set the result value and *signal* the monitor in order to trigger a status check cycle
            log.trace "Received task start ACK > ${message.memberId} to ${handler}"
            monitor.schedule {
                handler.runningMember = message.memberId
            }
        }
        else {
            // -- set an error result for the handler
            log.trace "Cluster member no more available: ${message.memberId} -- Notify error for taskId: ${message.taskId}"
            monitor.schedule {
                handler.result = HzCmdNotify.error( session.getUniqueId(), message.taskId )
            }
        }



//        if( handler.runningMember == null ) {
//            handler.status = TaskHandler.Status.COMPLETED
//        }
//        else {
//            log.trace "Ignoring handler: $handler -- this is supposed to be managed by 'memberRemoved' method"
//        }

    }

    @Override
    void itemRemoved(ItemEvent item) {
        /* ignored */
    }


    /**
     * Invoked by Hazelcast when a new cluster member is added
     *
     * @param membershipEvent
     */
    @Override
    void memberAdded(MembershipEvent e) {
        log.debug "Cluster member added: ${e.member}"

        HzNodeInfo node = allMembers.get( e.member.getUuid() )
        if( !node ) {
            log.debug "Oops .. no cluster node info for added member: ${e.member}"
            return
        }

        // -- increase the monitor capacity by the number of slots provided by the member added
        capacityInc(node.slots)

    }

    /**
     * Invoked by Hazelcast when a cluster member is removed
     *
     * @param membershipEvent
     */
    @Override
    void memberRemoved(MembershipEvent event) {
        log.debug "Cluster member removed > ${event.member} -- ${event.member?.getUuid()}"

        final memberId = event.member?.getUuid()
        final node = allMembers.get(memberId)
        if( !node ) {
            log.debug "Oops .. no node info for removed cluster: ${event.member}"
            return
        }

        // -- reduce the monitor capacity by the number of slots provided by the member removed
        capacityDec(node.slots)
        allMembers.remove(memberId)

        /*
         * !! reschedule all tasks that where running in a stopped node
         */
        def queue = new LinkedList<>(monitor.getPollingQueue()) as List<HzTaskHandler>
        def deadTasks = queue.findAll { HzTaskHandler handler ->
            log.trace "Checking status for handler: $handler"
            handler.runningMember == memberId && handler.status != TaskHandler.Status.COMPLETED
        }

        if( !deadTasks ) {
            log.trace "No dead tasks to be re-scheduled"
            return
        }

        log.debug "!! Some tasks were running on stopped node ${event.member} -- ${deadTasks}"
        monitor.schedule {
            for( HzTaskHandler handler : deadTasks ) {
                handler.result = HzCmdNotify.error( session.getUniqueId(), handler.task.id )
            }
        }

    }


    /**
     * @param member A {@code Member} instance
     * @return {@code true} when this member is in the map #allMembers, otherwise {@code false}
     */
    protected boolean isAvailable( Member member ) {
        isAvailable(member.getUuid())
    }

    /**
     * @param memberId A {@code Member} instance unique ID
     * @return {@code true} when this member is in the map #allMembers, otherwise {@code false}
     */
    protected boolean isAvailable( String memberId ) {
        assert memberId
        def result = allMembers.containsKey( memberId )
        return result
    }

    @Override
    void memberAttributeChanged(MemberAttributeEvent memberAttributeEvent) {
        /* ignore this */
    }


}
