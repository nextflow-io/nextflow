/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

import java.nio.channels.ClosedByInterruptException

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ProcessException
import nextflow.processor.TaskBean
import nextflow.processor.TaskId
import nextflow.processor.TaskRun
import nextflow.scheduler.Protocol.TaskResources
import nextflow.util.KryoHelper
import org.apache.ignite.compute.ComputeJob
import org.apache.ignite.lang.IgniteCallable
/**
 * Models a task executed remotely in a Ignite cluster node
 *
 * @param < T > The type of the value returned by the {@code #call} method
 */
@Slf4j
@CompileStatic
abstract class IgBaseTask<T> implements IgniteCallable<T>, ComputeJob {

    /**
     * Requested computational resources
     */
    TaskResources resources

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
    private TaskId taskId

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
        this.resources = new TaskResources(task)
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
     * Invoke the task execution. It calls the following methods in this sequence: {@code stage}, {@code execute0} and {@code unstage}
     *
     * @return The {@code execute0} result value
     * @throws nextflow.exception.ProcessException
     */
    @Override
    final T call() throws Exception {
        try {
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
        catch( InterruptedException | ClosedByInterruptException e ) {
            throw e
        }
        catch( Throwable e ) {
            log.error("Cannot execute task > ${bean?.name}", e)
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

    TaskId getTaskId() { taskId }

    @Override
    boolean equals( Object obj ) {
        if( obj.class != this.class ) return false
        final that = (IgBaseTask)obj
        return this.taskId == that.taskId
    }

    @Override
    int hashCode() {
        taskId.hashCode()
    }

    @Override
    String toString() {
        "${getClass().simpleName}[taskId=${taskId}]"
    }

}