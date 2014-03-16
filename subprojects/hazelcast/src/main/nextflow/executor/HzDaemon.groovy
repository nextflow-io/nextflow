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
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.ReentrantLock

import com.hazelcast.config.Config
import com.hazelcast.config.FileSystemXmlConfig
import com.hazelcast.config.UrlXmlConfig
import com.hazelcast.core.Hazelcast
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IQueue
import com.hazelcast.core.ITopic
import com.hazelcast.core.MembershipEvent
import com.hazelcast.core.MembershipListener
import com.hazelcast.instance.MemberImpl
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.daemon.DaemonLauncher
import nextflow.util.ConfigHelper

/**
 * Run the Hazelcast daemon used to process user processes
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzDaemon implements HzConst, DaemonLauncher, MembershipListener {

    private HazelcastInstance hazelcast

    private int capacity

    private ReentrantLock mutex

    private Condition notFull

    private volatile int count

    private IQueue<HzCmdCall> tasksQueue

    private ITopic<HzCmdResult> resultsTopic

    private ITopic<HzCmdStart> startsTopic

    private ExecutorService executor

    private Map<UUID, HzRemoteSession> allSessions

    /**
     * Poor man singleton. The current daemon instance
     */
    static HzDaemon current

    private Map config

    HzDaemon() {  }

    @Override
    void launch(Map properties) {
        log.info "Configuring Hazelcast cluster daemon"
        this.config = properties

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

    protected getHzProperty( String name, defValue = null ) {
        ConfigHelper.getConfigProperty(config, 'hazelcast', name, defValue)
    }

    @Memoized
    protected getHzPort () {
        getHzProperty('port') as Integer
    }

    @Memoized
    protected String getHzGroupName() {
        getHzProperty('group', DEFAULT_GROUP_NAME)
    }

    protected int getHzSlots() {
        getHzProperty('slots', Runtime.getRuntime().availableProcessors()) as Integer
    }

    protected List<String> getHzInterfaces() {
        def result = []
        def interfaceNames = getHzProperty('interface')?.toString()?.split(',') as List<String>
        interfaceNames?.each {
            if( it.contains('.') )
                result << it // it is supposed to be an interface IP address, add it to the list

            else
                result.addAll( findInterfaceAddressesByName(it) )
        }
        return result
    }

    protected List<String> findInterfaceAddressesByName( String name ) {
        assert name
        def result = []
        NetworkInterface interfaces = NetworkInterface.getNetworkInterfaces().find { NetworkInterface net -> net.name == name || net.displayName == name }
        interfaces?.getInetAddresses()?.each { InetAddress addr ->
            if( addr instanceof Inet4Address )
                result << addr.getHostAddress()
        }
        return result
    }

    @Memoized
    protected Config getConfigObj () {
        def result = new Config()
        def groupName = getHzGroupName()
        log.debug "Hazelcast config > group name: $groupName"
        result.getGroupConfig().setName(groupName)

        def network = result.getNetworkConfig()
        if( getHzPort() ) {
            log.debug "Hazelcast config > port: ${getHzPort()}"
            network.setPort( getHzPort() )
        }

        def interfaces = getHzInterfaces()
        if( interfaces ) {
            network.getInterfaces().setEnabled(true)
            interfaces.each {
                log.debug "Hazelcast config > setting network interface address: $it"
                network.getInterfaces().addInterface(it)
            }
        }

        def join = getHzProperty('join') as String ?: 'multicast'

        if( join == 'multicast') {
            log.debug "Hazelcast config > setup UDP/multicast"
            network.getJoin().getMulticastConfig().setEnabled(true)
            network.getJoin().getTcpIpConfig().setEnabled(false)
        }
        else {
            def members = join ? join.split(',') as List : []
            log.debug "Hazelcast config > setup with TCP joining members: $members"
            network.getJoin().getMulticastConfig().setEnabled(false)
            network.getJoin().getTcpIpConfig().setEnabled(true)

            if( members ) {
                network.getJoin().getTcpIpConfig().setMembers( members )
            }
        }

        log.debug("Hazelcast config > obj: ${result.toString()}")
        return result
    }

    /**
     * Initialize the daemon subsystem
     */
    protected initialize() {
        System.setProperty('hazelcast.logging.type','slf4j')
        System.setProperty('hazelcast.system.log.enabled','true')
        // -- see https://groups.google.com/d/msg/hazelcast/P-gu4em9WNk/b1uovn-k8rYJ
        System.setProperty('hazelcast.socket.bind.any', 'false')

        // -- define the number of jobs this manage in parallel
        capacity = getHzSlots()
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

        // -- set up distributed data structures
        tasksQueue = hazelcast.getQueue(TASK_SUBMITS_QUEUE)
        resultsTopic = hazelcast.getTopic(TASK_RESULTS_TOPIC)
        startsTopic = hazelcast.getTopic(TASK_STARTS_TOPIC)
        allSessions = hazelcast.getMap(SESSION_MAP)
        hazelcast.getCluster().addMembershipListener(this)

        // -- report info
        def cluster = hazelcast.getCluster()
        final local = cluster.getLocalMember() as MemberImpl
        final address = "${local.getAddress().getHost()}:${local.getAddress().getPort()}"
        log.info "Member [${address}] joined cluster '${cfg.getGroupConfig().getName()}' -- Total nodes: ${cluster.getMembers().size()}"

        // -- set the current instance
        current = this

        // -- shutdown on JVM exit
        Runtime.getRuntime().addShutdownHook {
            log.info "Shutting down node"
            try {
                hazelcast?.shutdown()
            }
            catch( Exception e ) {
                log.debug "Error during node shutdown", e
            }
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
        startsTopic.publish(new HzCmdStart(task))

        // -- execute the task
        executor.submit( {

            def result = null
            Throwable error = null
            try {
                result = task.call()
            }
            catch( Throwable ex ) {
                log.debug("!Error running task: $task", ex)
                error = ex
            }
            finally {
                publishCmdResult(task,result,error)
                mutex.withLock {
                    count--;
                    notFull.signal()
                }
            }


        } as Runnable )

    }

    /**
     * Publish the command result
     *
     * @param cmd The {@code HzCmdCall} that originates the result
     * @param result The actual command result
     * @param error The thrown exception eventually
     */
    private void publishCmdResult( HzCmdCall cmd, result, Throwable error ) {
        log.trace "Publishing command result: $cmd"

        try {
            final obj = new HzCmdResult(cmd,result,error)
            resultsTopic.publish(obj)

        }
        catch( Throwable throwable ) {
            log.debug("!Unable to publish command result: $cmd -- Publishing error cause instead: ${throwable.getMessage()?:throwable}")
            resultsTopic.publish( new HzCmdResult(cmd,result,throwable) )
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

    def String getHostNameAndAddress() {
        def host = InetAddress.getLocalHost()
        return "${host.getHostName()} [${host.getHostAddress()}]"
    }
}



