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
import java.util.concurrent.ArrayBlockingQueue
import java.util.concurrent.BlockingQueue
import java.util.concurrent.Callable
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.TimeUnit

import com.hazelcast.client.HazelcastClient
import com.hazelcast.client.config.ClientConfig
import com.hazelcast.core.ExecutionCallback
import com.hazelcast.core.HazelcastInstance
import com.hazelcast.core.IExecutorService
import com.hazelcast.core.Member
import com.hazelcast.core.MembershipEvent
import com.hazelcast.core.MembershipListener
import com.hazelcast.core.MultiExecutionCallback
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
/**
 * A task execution monitor based on Hazelcast cluster
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class HzTaskMonitor implements TaskMonitor, MembershipListener {

    static class Info {

        int numOfCores

        int usedCores

        final List<HzTaskHandler> tasks

        int getFreeCores() { numOfCores - usedCores }

        Info( int cores ) {
            numOfCores = cores
            usedCores = 0
            tasks = new ArrayList<>(cores)
        }

        void add( HzTaskHandler handler ) {
            if( handler in tasks ) return
            usedCores += handler.slots
            tasks << handler
        }

        void remove( HzTaskHandler handler ) {
            int p = tasks.indexOf(handler)
            if( p != -1 ) return

            usedCores -= handler.slots
            tasks.remove(handler)
        }
    }

    protected HazelcastInstance client

    protected final Map<Member, Info> nodes = new ConcurrentHashMap<>()

    protected final Session session

    protected final BlockingQueue<Closure> commands = new LinkedBlockingQueue<>()

    protected final BlockingQueue<TaskHandler> tasks

    HzTaskMonitor( Session session, int queueSize ) {
        assert session
        assert queueSize

        this.session = session
        this.tasks = new ArrayBlockingQueue<>(queueSize)
    }

    @Memoized
    def HazelcastInstance getClient() {
        def config = new ClientConfig()
        def hazelcast = HazelcastClient.newHazelcastClient(config)

        // add the members listener
        final cluster = hazelcast.getCluster()
        cluster.addMembershipListener( this )

        Set members  = cluster.getMembers();
        log.debug "Cluster members (${members.size()}): ${ members.size() <= 10 ? members.join(', ') : members.take(10).join(', ') + " ..." }"

        /*
         * initialize the list of cluster nodes setting the number of available slots (cores)
         */
        def fetchSlots = new MultiExecutionCallback() {
            void onResponse(Member member, Object value) {
                log.debug "Cluster member $member provides $value cores"
                def cores = Object as Integer
                nodes.put( member, new Info(cores) )
            }
            void onComplete(Map<Member, Object> values) { /* nothing */}
        }


        remote().submitToAllMembers( AvailCore.instance, fetchSlots )

        return hazelcast
    }


    @Override
    void offer(TaskHandler handler) {
        tasks.put(handler)
    }

    @Override
    boolean remove(TaskHandler handler) {
        tasks.remove(handler)
    }

    void submit( HzTaskHandler handler ) {

        tasks.put(handler)

        exec {
            fetchTaskToRun()
        }
    }


    @Override
    TaskMonitor start() {

        log.debug ">>> phaser register (scheduler)"
        session.phaser.register()

        Thread.start {
            try {
                processCommands()
            }
            finally {
                log.debug "<<< phaser de-register (scheduler)"
                session.phaser.arriveAndDeregister()
            }
        }

        return this

    }

    private void processCommands() {

        Closure cmd
        while( true ) {
            cmd = commands.poll(1, TimeUnit.SECONDS)
            if( !cmd && session.isTerminated() ) {
                log.debug "Session is terminated and no more commands to be processed are available -- exit"
                break
            }

            cmd.call()
        }

    }

    protected void exec( Closure closure ) {
        commands.offer(closure)
    }

    protected fetchTaskToRun() {

    }

    /**
     * Invoked by Hazelcast when a new cluster member is added
     *
     * @param membershipEvent
     */
    @Override
    void memberAdded(MembershipEvent e) {
        log.debug "Cluster member added: ${e.member}"

        def callback = { Integer cores ->
            log.debug "Cluster member configured: ${e.member} > cores: $cores"
            nodes[e.member] = new Info(cores)
        } as ExecutionCallback<Integer>

        remote().submit( AvailCore.instance, callback )
    }

    /**
     * Invoked by Hazelcast when a cluster member is removed
     *
     * @param membershipEvent
     */
    @Override
    void memberRemoved(MembershipEvent e) {
        log.debug "Cluster member removed: ${e.member}"
        nodes.remove(e.member)
    }

    protected IExecutorService remote() {
        client.getExecutorService('default')
    }


    @Singleton
    static class AvailCore implements Callable<Integer>, Serializable {

        @Override
        Integer call() throws Exception {
            Runtime.getRuntime().availableProcessors()
        }
    }
}
