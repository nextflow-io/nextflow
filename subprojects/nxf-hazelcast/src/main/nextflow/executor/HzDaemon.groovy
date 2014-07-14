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
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import com.hazelcast.config.Config
import com.hazelcast.config.FileSystemXmlConfig
import com.hazelcast.config.UrlXmlConfig
import com.hazelcast.core.EntryEvent
import com.hazelcast.core.EntryListener
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IMap
import com.hazelcast.core.IQueue
import com.hazelcast.core.Member
import com.hazelcast.core.MemberAttributeEvent
import com.hazelcast.core.MembershipEvent
import com.hazelcast.core.MembershipListener
import com.hazelcast.instance.MemberImpl
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.daemon.DaemonLauncher
import nextflow.util.DaemonConfig
import nextflow.util.Duration
import nextflow.util.RemoteSession
import org.apache.commons.lang.StringUtils
/**
 * Run the Hazelcast daemon used to process user processes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ServiceName('hazelcast')
class HzDaemon implements HzConst, DaemonLauncher, MembershipListener {

    static HEARTBEAT_DURATION = Duration.of('5s')

    private HazelcastInstance hazelcast

    private int capacity

    private ReentrantLock mutex

    private Condition notFull

    private volatile int count

    private IQueue<HzCmdCall> tasksQueue

    private ExecutorService executor

    private IMap<UUID, RemoteSession> allSessions

    private IMap<HzTaskKey, HzCmdStatus> allTasks

    private Map<String,HzNodeInfo> allNodes

    private Map<HzTaskKey, HzCmdStatus> runningTasks = [:]

    private final Lock tasksLock = new ReentrantLock()

    /**
     * Poor man singleton. The current daemon instance
     */
    static HzDaemon current

    private Member thisMember

    private DaemonConfig config

    HzDaemon() {  }

    @Override
    void launch(Map properties) {
        log.info "Configuring Hazelcast cluster daemon"
        config = new DaemonConfig('hazelcast', properties, System.getenv() )

        // -- initialize
        initialize()

        // -- wait for tasks to be processed
        processTasks()
    }

    /**
     * Read the configuration from a 'hazelcast.xml' file when available.
     * it looks for it in the following order:
     * 1) the URL/path specified by a the system property {@code hazelcast.config}
     * 2) the current folder
     * 3) the nextflow home folder
     *
     * @return The {@code Config} object when the file is available or {@code null} otherwise
     */
    protected Config tryNativeConfig() {

        final localConf = new File('hazelcast.xml')
        final homeConf = new File(Const.APP_HOME_DIR, 'hazelcast.xml')

        // -- check for an external config file
        Config cfg = null
        String confFile = System.getProperty('hazelcast.config')

        if( !confFile && localConf.exists() )
            confFile = localConf.toString()

        if( !confFile && homeConf.exists() )
            confFile = homeConf.toString()


        if( confFile ) {
            try {
                log.debug "Hazelcast config > loading native config: $confFile"
                cfg = confFile.contains('://') ? new UrlXmlConfig(confFile) : new FileSystemXmlConfig(confFile)
            }
            catch( MalformedURLException e ) {
                log.warn "Hazelcast config > Unable to read native configuration: $confFile -- Make sure you entered a valid URL"
            }
            catch( FileNotFoundException e ) {
                log.warn "Hazelcast config > Unable to read native configuration: $confFile -- Make sure the file you specified exist"
            }
            catch( Exception e ) {
                log.warn "Hazelcast config > Unable to read native configuration: $confFile -- Cause: ${e.getMessage()}"
            }
        }

        return cfg
    }


    /**
     * @return The daemon port
     */
    protected getPort() {
        config.getAttribute('port') as Integer
    }

    /**
     * @return The cluster group name
     */
    protected String getGroupName() {
        config.getAttribute('group', DEFAULT_GROUP_NAME)
    }

    /**
     * @return The slot provided by the cluster node i.e. the number of processor
     */
    protected int getSlots() {
        config.getAttribute('slots', Runtime.getRuntime().availableProcessors()) as Integer
    }

    /**
     * @link http://hazelcast.org/docs/latest/manual/html-single/#system-property
     *
     * @return Minimum interval to consider a connection error as critical in milliseconds.
     */
    protected String getConnectionMonitorInterval() {
        config.getAttribute('connectionMonitorInterval', '300') as String
    }

    /**
     * @link http://hazelcast.org/docs/latest/manual/html-single/#system-property
     *
     * @return Maximum IO error count before disconnecting from a node
     */
    protected String getConnectionMonitorMaxFaults() {
        config.getAttribute('connectionMonitorMaxFaults', '5') as String
    }


    /**
     * @return Creates the hazelcast configuration object
     */
    @Memoized
    protected Config getConfigObj () {
        def result = new Config()
        def groupName = groupName
        log.debug "Hazelcast config > group name: $groupName"
        result.getGroupConfig().setName(groupName)

        def network = result.getNetworkConfig()
        if( getPort() ) {
            log.debug "Hazelcast config > port: ${port}"
            network.setPort(port)
        }

        def interfaces = config.getNetworkInterfaceAddresses()
        if( interfaces ) {
            network.getInterfaces().setEnabled(true)
            interfaces.each {
                log.debug "Hazelcast config > setting network interface address: $it"
                network.getInterfaces().addInterface(it)
            }
        }

        def join = config.getAttribute('join') as String ?: 'multicast'

        if( join == 'multicast') {
            log.debug "Hazelcast config > setup UDP/multicast"
            network.getJoin().getMulticastConfig().setEnabled(true)
            network.getJoin().getTcpIpConfig().setEnabled(false)
        }
        else {
            def members = join ? StringUtils.split(join, ", \n").collect { it.trim() } : []
            log.debug "Hazelcast config > setup with TCP joining members: $members"
            network.getJoin().getMulticastConfig().setEnabled(false)
            network.getJoin().getTcpIpConfig().setEnabled(true)

            if( members ) {
                network.getJoin().getTcpIpConfig().setMembers( members )
            }
        }

        // -- config data structures
        //result.getMapConfig(MEMBERS_MAP)

        log.debug("Hazelcast config > obj: ${result.toString()}")
        return result
    }

    /**
     * Initialize the daemon subsystem
     */
    protected initialize() {
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')
        System.setProperty('hazelcast.memcache.enabled','false')
        System.setProperty('hazelcast.rest.enabled', 'false')
        System.setProperty('hazelcast.version.check.enabled', 'false')
        // -- see https://groups.google.com/d/msg/hazelcast/P-gu4em9WNk/b1uovn-k8rYJ
        System.setProperty('hazelcast.socket.bind.any', 'false')

        // -- connection properties
        final interval = getConnectionMonitorInterval()
        log.debug "Hazelcast config > connectionMonitorInterval: $interval millis"
        System.setProperty('hazelcast.connection.monitor.interval', interval)

        final maxFaults = getConnectionMonitorMaxFaults()
        log.debug "Hazelcast config > connectionMonitorMaxFaults: $maxFaults"
        System.setProperty('hazelcast.connection.monitor.max.faults', maxFaults)

        // -- define the number of jobs this manage in parallel
        capacity = getSlots()
        log.debug "Hazelcast config > slots: $capacity"
        if( capacity <= 0 )
            throw new IllegalArgumentException("Daemon 'slots' parameter cannot be less than 1")

        // -- configure the hazelcast daemon
        def cfg = tryNativeConfig()
        if( !cfg ) {
            cfg = getConfigObj()
        }

        // -- add serializer and create a new instance
        HzSerializerConfig.registerAll(cfg.getSerializationConfig())
        hazelcast = Hazelcast.newHazelcastInstance(cfg)
        thisMember = hazelcast.getCluster().getLocalMember()

        // -- set up distributed data structures
        tasksQueue = hazelcast.getQueue(TASK_SUBMITS_QUEUE)
        allNodes = hazelcast.getMap(ALL_NODES_MAP)
        allTasks = hazelcast.getMap(ALL_TASKS_MAP)
        allSessions = hazelcast.getMap(SESSIONS_MAP)
        allSessions.addEntryListener( new SessionEntryListener(), false )
        hazelcast.getCluster().addMembershipListener(this)

        // -- report info
        def cluster = hazelcast.getCluster()
        final local = thisMember as MemberImpl
        final address = "${local.getAddress().getHost()}:${local.getAddress().getPort()}"
        log.info "Member [${address}] joined cluster '${cfg.getGroupConfig().getName()}' -- Total nodes: ${cluster.getMembers().size()}"

        // -- set the current node info
        allNodes.put( local.getUuid(), new HzNodeInfo(capacity)  )

        // -- set the current instance
        current = this

        registerHeartbeat()
        // shutdown on termination
        Runtime.getRuntime().addShutdownHook { shutdown() }
    }

    /**
     * launch a thread that update periodically the current running tasks
     */
    final protected void registerHeartbeat() {
        final MAX_DELTA = HEARTBEAT_DURATION.toMillis()

        def heartbeat = {

            while( true ) {
                sleep(MAX_DELTA)
                if( !hazelcast?.getLifecycleService()?.isRunning() ) break

                final now = System.currentTimeMillis()
                runningTasks.each { key, value ->
                    if( now-value.lastUpdate < MAX_DELTA )
                        return

                    log.debug "% Update status for task: $value"
                    def newValue = value.copyWith( clock: now )

                    tasksLock.with {
                        runningTasks.put(key, newValue)
                        def replaced = allTasks.replace(key, value, newValue)
                        if( !replaced ) {
                            log.warn "% Failed to update task status: $newValue"
                        }
                    }

                }
            }
            log.debug "% Heartbeat thread terminated"

        } as Runnable


        def thread = new Thread(heartbeat)
        thread.setName('HzDaemon heartbeat')
        thread.setDaemon(true)
        thread.run()
    }

    /**
     * Shutdown the Hazelcast instance
     */
    final protected void shutdown() {
        if( !hazelcast ) return
        log.info "Shutting down node"
        try {
            hazelcast.shutdown()
            hazelcast = null
        }
        catch( Throwable e ) {
            // do not care
        }
    }

    /**
     * Get a class loader for a remote script
     *
     * @param sessionId The client identifier
     * @return The class loader associated to the remote client
     */
    @Memoized
    static ClassLoader getClassLoaderFor(UUID sessionId) {
        if( !current )
            throw new IllegalStateException('Missing daemon instance')

        def session = current.allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing remote session for ID: $sessionId")

        def loader = new GroovyClassLoader()
        session.classpath.each { URL url ->
            log.debug "Adding classpath url: $url"
            loader.addURL(url)
        }

        return loader
    }

    /**
     * Takes a task to be processed from the distributed queue an execute it
     */
    protected void processTasks() {

        // the shared queue
        mutex = new ReentrantLock()
        notFull = mutex.newCondition()
        executor = Executors.newCachedThreadPool()

        Thread.start {
            count = 0

            // make sure to do not run more tasks than that the actual capacity
            while( true ) {
                mutex.withLock {
                    while( count >= capacity ) { notFull.await() }
                    count++
                }

                execute(tasksQueue.take())
            }

        }
    }

    /**
     * Executes the received command and send back the result when finished
     *
     * @param task The command task to be executed
     */
    protected void execute( HzCmdCall task ) {
        log.trace "Executing task command > $task"

        // -- notify the execution of this task
        notifyStart(task)

        // -- execute the task
        executor.submit( {

            def result = null
            Throwable error = null
            try {
                result = task.call()
            }
            catch( Throwable ex ) {
                log.debug("Error running task: $task", ex)
                error = ex
            }
            finally {
                notifyResult(task,result,error)
                mutex.withLock {
                    count--;
                    notFull.signal()
                }
            }


        } as Runnable )

    }

    /**
     * Add the command start notification in the event queue
     * @param cmd
     */
    private void notifyStart( HzCmdCall cmd ) {

        // add this task id in the map of all processed tasks
        final stateObj = HzCmdStatus.start(cmd, thisMember.getUuid())
        final key = new HzTaskKey( cmd.sessionId, cmd.taskId )

        tasksLock.withLock {
            // map of current running tasks
            runningTasks[ key ] = stateObj
            // add to the shared obj
            allTasks[ key ] = stateObj
        }

    }

    /**
     * Publish the command result
     *
     * @param cmd The {@code HzCmdCall} that originates the result
     * @param result The actual command result
     * @param error The thrown exception eventually
     */
    private void notifyResult( HzCmdCall cmd, result, Throwable error ) {
        log.trace "Publishing command result: $cmd"

        // add this task id in the map of all processed tasks
        final key = new HzTaskKey( cmd.sessionId, cmd.taskId )

        tasksLock.withLock {
            allTasks[ key ] = HzCmdStatus.result(cmd, result, error)
            runningTasks.remove( key )
        }

    }

    @Override
    void memberAdded(MembershipEvent event) {
        log.info "Nextflow cluster member added: ${event.member}"
    }

    @Override
    void memberRemoved(MembershipEvent event) {
        log.info "Nextflow cluster member remove: ${event.member}"
    }

    @Override
    void memberAttributeChanged(MemberAttributeEvent memberAttributeEvent) {
        /* ignore this */
    }


    @Slf4j
    class SessionEntryListener implements EntryListener<UUID,RemoteSession> {

        /**
         * When a session entry is removed (i.e. when a client disconnect)
         * remove all the eventually remaining task id
         * @param event
         */
        @Override
        void entryRemoved(EntryEvent<UUID, RemoteSession> event) {
            final sessionId = event.getKey()
            log.debug "Removing session id: ${sessionId}"
            def list = []
            HzDaemon.this.allTasks.keySet().each { HzTaskKey key ->
                if( key.sessionId == sessionId ) {
                    list << key.taskId
                    HzDaemon.this.allTasks.remove(key)
                }
            }
            if( list )
                log.debug "Evicted ${list.size()} tasks from session: $sessionId -- $list"
        }

        /** do not used */
        @Override
        void entryAdded(EntryEvent<UUID, RemoteSession> event) { }

        /** do not used */
        @Override
        void entryUpdated(EntryEvent<UUID, RemoteSession> event) { }

        /** do not used */
        @Override
        void entryEvicted(EntryEvent<UUID, RemoteSession> event) { }
    }
}



