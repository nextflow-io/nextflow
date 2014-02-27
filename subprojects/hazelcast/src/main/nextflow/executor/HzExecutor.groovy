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
import java.nio.file.Path
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
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.script.ScriptType
import nextflow.util.Duration
/**
 * Creates an executor that submits process to an Hazelcast cluster for execution
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HzExecutor extends AbstractExecutor  {

    @PackageScope
    HzConnector connector

    /**
     * Initialize the executor by getting a reference to the Hazelcast connector
     */
    def void init() {
        super.init()
        connector = HzConnector.create(taskMonitor)
    }

    /**
     * Creates the task monitor for this executor
     * @return An instance of {@code TaskMonitor}
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, Duration.of('1s'))
    }


    /**
     *  Creates an handler for the specified task
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {

        if( task.type == ScriptType.GROOVY ) {
            throw new UnsupportedOperationException("Native tasks are not supported")
        }

        /*
         * otherwise as a bash script
         */
        final bash = new BashWrapperBuilder(task)
        bash.environment = task.processor.getProcessEnvironment()
        bash.environment.putAll( task.getInputEnvironment() )

        // staging input files
        bash.stagingScript = {
            final files = task.getInputFiles()
            final staging = stagingFilesScript(files)
            return staging
        }

        // unstage script
        bash.unstagingScript = {
            return unstageOutputFilesScript(task)
        }

        // create the wrapper script
        bash.build()

        new HzTaskHandler(task, taskConfig, this)
    }


    TaskPollingMonitor getTaskMonitor() {
        (TaskPollingMonitor)super.getTaskMonitor()
    }


    /**
     * The sender ID is the current session ID
     */
    @PackageScope
    UUID getSenderId() {
        return session.uniqueId
    }

}

/**
 * Creates a connector for the Hazelcast cluster
 *
 */
@Slf4j
class HzConnector implements HzConst, MembershipListener, MessageListener<HzBashResult> {

    private HazelcastInstance hazelcast

    private ITopic<HzBashResult> resultsTopic

    private IQueue<HzBashCmd> executorsQueue

    @PackageScope
    private final Map<Member,Integer> slotsFor = [:]


    private TaskPollingMonitor monitor

    /** The nextflow session object */
    private Session session

    private Map<UUID,RemoteSession> allSessions


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
    Queue<HzBashCmd> getExecutorsQueue() {
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
        allSessions.put(session.uniqueId, new RemoteSession(session)  )

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
     * @param message The {@code HzBashResult} represent the result of a task executed in the Hazelcast cluster
     */
    @Override
    void onMessage(Message<HzBashResult> message) {
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

/**
 * A task handler for Hazelcast cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzTaskHandler extends TaskHandler {

    private HzExecutor executor

    private final Path exitFile

    private final Path wrapperFile

    private final Path outputFile

    @PackageScope
    volatile HzBashResult result

    protected HzTaskHandler(TaskRun task, TaskConfig taskConfig, HzExecutor executor) {
        super(task, taskConfig)
        this.executor = executor
        this.exitFile = task.getCmdExitFile()
        this.outputFile = task.getCmdOutputFile()
        this.wrapperFile = task.getCmdWrapperFile()
    }

    @Override
    void submit() {
        final List cmd = new ArrayList(taskConfig.shell ?: 'bash' as List ) << wrapperFile.getName()

        // submit to an hazelcast node for execution
        def command = new HzBashCmd( executor.senderId, task.id, task.workDirectory, cmd )
        executor.connector.executorsQueue.add( command )

        // mark as submitted -- transition to STARTED has to be managed by the scheduler
        status = TaskHandler.Status.SUBMITTED
    }

    @Override
    boolean checkIfRunning() {
        if( isSubmitted() ) {
            log.trace "Task ${task.name} > RUNNING"
            status = TaskHandler.Status.RUNNING
            return true
        }

        return false
    }

    @Override
    boolean checkIfCompleted() {

        if( isRunning() && result != null ) {
            // TODO Add a check on the result file existence due to NFS latency
            status = TaskHandler.Status.COMPLETED

            // -- set the task exit code (only when it is a scriptlet task)
            if( result.isScriptlet() && result.value instanceof Integer )
                task.exitStatus = result.value as Integer

            // -- the task output depend by the kind of the task executed
            if( result.isScriptlet() )
                task.stdout = outputFile
            else
                task.stdout = result.value

            // -- set the task result error (if any)
            task.error = result.error

            log.trace "Task ${task.name} > DONE (${task.exitStatus})"
            return true
        }

        return false
    }

    @Override
    void kill() {
        // not implemented
        // see also https://groups.google.com/d/msg/hazelcast/-SVy4k1-QrA/rlUee3i4POAJ
    }

}

