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
import com.hazelcast.core.EntryEvent
import com.hazelcast.core.EntryListener
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IExecutorService
import com.hazelcast.core.IMap
import com.hazelcast.core.IQueue
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
class HzConnector implements HzConst, MembershipListener {

    private HazelcastInstance hazelcast

    private IQueue<HzCmdCall> submitQueue

    /** The nextflow session object */
    private Session session

    private IMap<UUID,HzRemoteSession> allSessions

    private IMap<String,HzNodeInfo> allNodes

    private IMap<HzTaskKey, HzCmdStatus> allTasks

    private TaskPollingMonitor monitor

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

    @PackageScope
    IExecutorService remote() {
        hazelcast.getExecutorService(EXEC_SERVICE)
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

        allTasks = hazelcast.getMap(ALL_TASKS_MAP)
        allTasks.addEntryListener( new TaskEntryListener(), true )

        // publish the session on the cluster
        allSessions = hazelcast.getMap(SESSIONS_MAP)
        allSessions.put( session.uniqueId, new HzRemoteSession(session) )

        // fetch the current number of slots for each member
        def tot = 0
        allNodes = hazelcast.getMap(ALL_NODES_MAP)
        allNodes.each { key, node -> tot += node.slots }
        capacityInc(tot)

        // add the members listener
        final cluster = hazelcast.getCluster()
        cluster.addMembershipListener( this )
        Set members  = cluster.getMembers();
        log.debug "Cluster members (${members.size()}): ${ members.size() <= 10 ? members.join(', ') : members.take(10).join(', ') + " ..." }"

        monitor.register( this.&pollingWatchdog )
        monitor.dispatcher.addCompleteListener(this.&removeCompleteTask)

        // register for the shutdown
        session.onShutdown {
            if( !hazelcast ) return
            log.info "Shutting down hazelcast connector"
            allSessions.remove(session.getUniqueId())
            allTasks.keySet().each { HzTaskKey key -> if(key.sessionId==session.getUniqueId()) allTasks.remove(key) }
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

    @PackageScope
    void setMember( Member member, HzNodeInfo info ) {
        assert info
        allNodes[member.getUuid()] = info
        capacityInc( info.slots )
    }


    /**
     * Invoked by Hazelcast when a new cluster member is added
     *
     * @param membershipEvent
     */
    @Override
    void memberAdded(MembershipEvent e) {
        log.debug "Cluster member added: ${e.member}"

        HzNodeInfo node = allNodes.get( e.member.getUuid() )
        if( node ) {
            // -- increase the monitor capacity by the number of slots provided by the member added
            capacityInc(node.slots)
        }
        else {
            log.debug "Oops .. no cluster node info for added member: ${e.member}"
            int sum = allNodes.values() *. slots .sum() as int
            int cnt = allNodes.size()
            int avg = cnt && sum ? sum/cnt : 1
            capacityInc(avg)
        }

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
        final node = allNodes.get(memberId)
        if( node ) {
            // -- reduce the monitor capacity by the number of slots provided by the member removed
            capacityDec(node.slots)
            allNodes.remove(memberId)
        }
        else {
            log.debug "Oops .. no node info for removed cluster: ${event.member}"
        }


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
                handler.result = HzCmdStatus.error( session.getUniqueId(), handler.task.id )
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
        def result = allNodes.containsKey( memberId )
        return result
    }

    @Override
    void memberAttributeChanged(MemberAttributeEvent memberAttributeEvent) {
        /* ignore this */
    }

    final protected HzTaskKey key( taskId ) {
        new HzTaskKey( session.getUniqueId(), taskId )
    }


    final protected void removeCompleteTask( TaskHandler handler ) {
        if( handler instanceof HzTaskHandler ) {
            log.debug "Removing from AllTasks taskId: ${handler.task.id}"
            allTasks.remove( key(handler.task.id) )
        }
    }

    /**
     * Watch the scheduled task status
     */
    final protected void pollingWatchdog( TaskPollingMonitor monitor ) {

        // find all still in submit status after more that X
        // if not in the queue --> bad
        // if it is in the queue and there is a free slot --> bad

        final now = System.currentTimeMillis()
        final interval = HzDaemon.HEARTBEAT_DURATION.toMillis() * 5

        final handlers = (List<HzTaskHandler>) monitor.findTaskHandlers { HzTaskHandler handler -> handler.isSubmitted() && now-handler.lastUpdate > interval }
        final pendingIds = handlers.collect { it.task.id }
        if( pendingIds ) {
            final queueIds = submitQueue.findAll { HzCmdCall it -> it.sessionId == session.uniqueId }.collect { HzCmdCall it -> it.taskId }
            for( def taskId : pendingIds ) {
                if( taskId in queueIds ) {
                    // the taskId in the queue of submitted tasks, it should check if all nodes are busy
                    // otherwise there's something bad
                    // TODO
                    log.debug "% Task is in the submit queue waiting to be processed > taskId: $taskId"
                }
                else {
                    log.debug "% Task in submit status but does not exist in the submit queue > taskId: $taskId"
                    def bad = handlers.find { HzTaskHandler it -> it.task.id == taskId }
                    if( bad ) {
                        bad.result = HzCmdStatus.error( session.getUniqueId(), taskId )
                    }
                    else {
                        log.debug "% Missing handler for submitted taskId: $taskId"
                    }
                }
            }
        }


        // find all in running state and last update more than X
        // --> kill and resubmit
        final runningTasks = allTasks.values().findAll {
            HzCmdStatus it -> it.sessionId == session.getUniqueId() && it.isStart() && now-it.lastUpdate > interval
        }

        if( runningTasks ) {
            for( HzCmdStatus cmd : runningTasks ) {
                def handler = (HzTaskHandler)monitor.getTaskHandlerById(cmd.taskId)
                if( handler ) {
                    log.debug "% Task in running status evicted > taskId: ${cmd.taskId}"
                    handler.result = HzCmdStatus.error( session.getUniqueId(), cmd.taskId )
                }
                else {
                    log.debug "% Missing handler for runnig taskId: ${cmd.taskId}"
                }
            }
        }

    }


    @Slf4j
    class TaskEntryListener implements EntryListener<HzTaskKey, HzCmdStatus> {

        @Override
        void entryAdded(EntryEvent<HzTaskKey, HzCmdStatus> event) {
            if( event.key.sessionId == session.uniqueId ) {
                log.trace "Received ADD task event: $event"
                handleEntry( event.getValue() )
            }
            else {
                log.trace "Discarding ADD task event for other client: ${event}"
            }
        }

        @Override
        void entryUpdated(EntryEvent<HzTaskKey, HzCmdStatus> event) {
            if( event.key.sessionId == session.uniqueId ) {
                log.trace "Received UPDATE task event: $event"
                handleEntry( event.getValue() )
            }
            else {
                log.trace "Discarding UPDATE task event for other client: ${event}"
            }
        }

        @Override
        void entryRemoved(EntryEvent<HzTaskKey, HzCmdStatus> event) {
            /* ignore this */
        }


        @Override
        void entryEvicted(EntryEvent<HzTaskKey, HzCmdStatus> event) {
            log.debug "AllTasks evicted entry  > ${event.getValue()}"
        }

        private handleEntry ( HzCmdStatus event ) {

            if( !event ) {
                log.trace "Oops .. Empty task event -- ignore it"
                return
            }

            /*
             * find out the task handler for the task id of the completed task
             */
            def handler = (HzTaskHandler)monitor.getTaskHandlerById(event.taskId)
            if( !handler ) {
                log.warn "Unknow task with ID: ${event.taskId}"
                return
            }

            /*
             * when it is a task completion notification,
             * update the handler with that object
             */
            if( event.isComplete() ) {
                log.trace "Task complete > ${event}"
                monitor.schedule {
                    handler.result = event
                }
                return // exit here
            }

            /*
             * Otherwise this event notifies that the task has started
             * Make sure the cluster member is still available
             */
            if( isAvailable(event.memberId) ) {
                // -- set the result value and *signal* the monitor in order to trigger a status check cycle
                log.trace "Received task start ACK > ${event.memberId} to ${handler}"
                monitor.schedule {
                    handler.runningMember = event.memberId
                }
            }
            else {
                // -- set an error result for the handler
                log.trace "Cluster member no more available: ${event.memberId} -- Notify error for taskId: ${event.taskId}"
                monitor.schedule {
                    handler.result = HzCmdStatus.error( session.getUniqueId(), event.taskId )
                }
            }
        }


    }

}
