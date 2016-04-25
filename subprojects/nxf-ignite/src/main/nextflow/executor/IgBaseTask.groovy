/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.CompileStatic
import nextflow.scheduler.JobComputeResources
import nextflow.daemon.IgGridFactory
import nextflow.exception.ProcessException
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.KryoHelper
import nextflow.util.RemoteSession
import org.apache.ignite.Ignite
import org.apache.ignite.IgniteCache
import org.apache.ignite.IgniteLogger
import org.apache.ignite.compute.ComputeJob
import org.apache.ignite.lang.IgniteCallable
import org.apache.ignite.resources.IgniteInstanceResource
import org.apache.ignite.resources.LoggerResource
/**
 * Models a task executed remotely in a Ignite cluster node
 *
 * @param < T > The type of the value returned by the {@code #call} method
 */
@CompileStatic
abstract class IgBaseTask<T> implements IgniteCallable<T>, ComputeJob {

    static final Map<UUID,GroovyClassLoader> classLoaderCache = new HashMap()

    @LoggerResource
    private transient IgniteLogger log

    @IgniteInstanceResource
    private transient Ignite grid

    /**
     * Requested computational resources
     */
    JobComputeResources resources

    /**
     * The client session identifier, it is required in order to access to
     * remote class-path
     */
    UUID sessionId

    /**
     * This field is used to transport the class attributes as a unique serialized byte array
     */
    private byte[] payload

    /**
     * Task unique identifier
     */
    private taskId

    /**
     * Holds the class attributes in this map. Note: is defined as 'transient' because
     * the map content is serialized as a byte[] and saved to the {@code payload} field
     */
    protected transient TaskBean bean

    /**
     * Initialize the grid gain task wrapper
     *
     * @param task The task instance to be executed
     * @param The session unique identifier
     */
    protected IgBaseTask( TaskRun task, UUID sessionId ) {
        this.sessionId = sessionId
        this.taskId = task.id
        this.bean = new TaskBean(task)
        this.payload = KryoHelper.serialize(bean)
        this.resources = new JobComputeResources(task)
    }

    /** ONLY FOR TESTING PURPOSE */
    protected IgBaseTask() {}

    /**
     * Template method invoked before the task is executed
     */
    protected void beforeExecute() {

    }

    /**
     * Template method invoke after the task is executed
     */
    protected void afterExecute() {

    }

    final protected void deserialize() {

        if( bean == null && payload ) {
            bean = (TaskBean)KryoHelper.deserialize(payload)
        }

    }

    /**
     * Notify that the task execution has started
     */
    private void notifyStart() {
        def master = grid.cluster().forAttribute(IgGridFactory.NODE_ROLE, IgGridFactory.ROLE_MASTER)
        grid.message(master).send(IgConnector.TOPIC_EVT_TASKS, taskId)
    }

    /**
     * Invoke the task execution. It calls the following methods in this sequence: {@code stage}, {@code execute0} and {@code unstage}
     *
     * @return The {@code execute0} result value
     * @throws nextflow.exception.ProcessException
     */
    @Override
    final T call() throws Exception {
        try {
            notifyStart()
            deserialize()

            /*
             * stage the input files in the working are`
             */
            beforeExecute()

            /*
             * execute the task
             */
            final T result = execute0()

            /*
             * copy back the result files to the shared area
             */
            afterExecute()

            // return the exit status eventually
            return result
        }
        catch( Exception e ) {
            log.error("Cannot execute task > $bean.name", e)
            throw new ProcessException(e)
        }

    }

    /**
     * Just a synonym for {@code #call}
     *
     * @return The value returned by the task execution
     */
    final Object execute() {
        call()
    }

    /**
     * The actual task executor code provided by the extending subclass
     *
     * @return The value returned by the task execution
     */
    protected abstract T execute0()

    /**
     * Lookup the {@link RemoteSession} object for the given session ID
     *
     * @param sessionId The remote session ID
     * @param grid The current {@link Ignite} instance
     * @return The associated {@link RemoteSession} object for the specified session ID
     * @throws IllegalStateException when no session is found for the specified session ID
     */
    protected RemoteSession getSessionFor( UUID sessionId ) {
        assert sessionId
        IgniteCache<UUID, RemoteSession> allSessions = grid.cache( IgGridFactory.SESSIONS_CACHE )

        if( !allSessions )
            throw new IllegalStateException('Missing session cache object')

        def session = allSessions.get(sessionId)
        if( !session )
            throw new IllegalStateException("Missing session object for id: $sessionId")

        return session
    }

    def getTaskId() { taskId }

    @Override
    String toString() {
        "${getClass().simpleName}[taskId: ${taskId}]"
    }

    /**
     * Create a {@link ClassLoader} object for the specified session ID
     *
     * @param sessionId
     * @param grid
     * @return
     */
    protected ClassLoader getClassLoaderFor( UUID sessionId ) {
        assert sessionId

        classLoaderCache.getOrCreate(sessionId) {

            final allSessions = (IgniteCache<UUID, RemoteSession>)grid.cache( IgGridFactory.SESSIONS_CACHE )
            if( !allSessions )
                throw new IllegalStateException('Missing session cache object')

            final session = allSessions.get(sessionId)
            if( !session )
                throw new IllegalStateException("Missing session object for id: $sessionId")

            final result = new GroovyClassLoader()
            session.classpath.each { Path file ->
                log.debug "Adding to classpath: $file"
                result.addClasspath(file.toAbsolutePath().toString())
            }

            return result
        }
    }
}