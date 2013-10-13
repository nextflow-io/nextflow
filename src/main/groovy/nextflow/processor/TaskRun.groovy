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

package nextflow.processor

import java.nio.file.Path

import groovy.transform.ToString
import nextflow.executor.ProcessHandler

/**
 * Models a task instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@ToString( includePackage = false, includeNames = true, includes = 'id,index,name,status,exitCode' )
class TaskRun<T extends ProcessHandler> {

    enum Status { NEW, STARTED, TERMINATED }

    /**
     * Task unique id
     */
    def id

    /**
     * TODO this have to be refactored
     * The jobid as provided by the underlying execution system (process, SGE jobid, dnanexus jobid, etc)
     */
    def jobId

    /**
     * Task index within its execution group
     */
    def index

    /**
     * Task name
     */
    def String name

    /*
     * The processor that creates this 'task'
     */
    TaskProcessor processor

    /**
     * The process handler which is executing the task
     */
    def T handler

    /**
     * Task current status
     */
    Status status

    /**
     * Holds the input value(s) for each task input parameter
     */
    Map<InParam,Object> inputs = [:]

    /**
     * Holds the output value(s) for each task output parameter
     */
    Map<OutParam,Object> outputs = [:]

    /**
     * Return the list of all input files staged as inputs by this task execution
     */
    @Lazy
    List<Path> stagedInputs = {
        stagedProvider.call().values().flatten() *. stagePath
    } ()

    /**
     * The *strategy* sued to retrieve the list of staged input files for this task.
     * This sort of *hack* is required since task processed by the {@code MergeTaskProcessor} maintain
     * their own input files list, and so the task will need to access that list, and not the one
     * hold by the task itself
     *
     * See MergeTaskProcessor
     */
    Closure<Map<FileInParam,List<FileHolder>>> stagedProvider = {
           (Map<FileInParam,List<FileHolder>>) getInputsByType(FileInParam)
    }


    def void setInput( InParam param, Object value = null ) {
        assert param

        inputs[param] = value

        // copy the value to the task 'input' attribute
        // it will be used to pipe it to the process stdin
        if( param instanceof StdInParam) {
            stdin = value
        }
    }

    def void setOutput( OutParam param, Object value = null ) {
        assert param

        outputs[param] = value
        if( param instanceof StdOutParam ) {
            stdout = value
        }
    }

    /**
     * The value to be piped to the process stdin
     */
    def stdin

    /**
     * The exit code returned by executing the task script
     */
    int exitCode = Integer.MAX_VALUE

    /**
     * The return std out
     */
    def stdout

    def String getStdout() {

        if( stdout instanceof Path ) {
            return stdout.text
        }
        else {
            stdout?.toString()
        }

    }


    /**
     * The directory used to run the task
     */
    Path workDirectory

    /*
     * The closure implementing this task
     */
    Closure code

    /**
     * The script executed by the system
     */
    def script

    def String getScript() {
        if( script instanceof Path ) {
            return script.text
        }
        else {
            return script?.toString()
        }
    }

    def void setStatus( Status status ) {
        assert status
        if ( this.status == status ) { return }

        this.status = status
        if( status == Status.NEW ) {
            processor.notifyTaskCreated(this)
        }
        if ( status == Status.STARTED ) {
            processor.notifyTaskStarted(this)
        }
        else if ( status == Status.TERMINATED ) {
            processor.notifyTaskCompleted(this)
        }
    }

    boolean isNew() { return status == Status.NEW }

    boolean isStarted() { return status == Status.STARTED }

    boolean isTerminated()  { return status == Status.TERMINATED  }

    /** The timestamp when the task has been submitted for execution */
    long submitTimeMillis

    /** The timestamp when the task has started the execution */
    long startedTimeMillis

    /** The timestamp when the task has completed the execution */
    long completedTimeMillis

    def <T extends InParam> Map<T,Object> getInputsByType( Class<T>... types ) {

        def result = [:]
        inputs.findAll() { types.contains(it.key.class) }.each { result << it }
        return result
    }

    def Map<OutParam,Object> getOutputsByType( Class... types ) {
        def result = [:]
        outputs.findAll() { types.contains(it.key.class) }.each { result << it }
        return result
    }



}

