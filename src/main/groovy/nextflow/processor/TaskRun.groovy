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
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.ToString
import groovy.util.logging.Slf4j

/**
 * Models a task instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ToString( includePackage = false, includeNames = true, includes = 'id,index,name,status,exitCode' )
class TaskRun<T extends TaskHandler> {

    static private final echoLock = new ReentrantLock(true)

    /**
     * Task unique id
     */
    def id

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

    /**
     * @return The task produced stdout result as string
     */
    def String getStdout() {

        if( stdout == null ) {
            stdout = getCmdOutputFile()
        }

        if( stdout instanceof Path ) {
            return stdout.text
        }
        else {
            return stdout?.toString()
        }
    }

    /**
     * Print to the current system console the task produced stdout
     */
    void echoStdout() {

        def out = stdout ?: getCmdOutputFile()

        // print the stdout
        if( out instanceof Path ) {
            if( !out.exists() ) {
                log.debug "Echo file does not exist: ${out}"
                return
            }

            echoLock.lock()
            try {
                out.withReader {  System.out << it }
            }
            finally {
                echoLock.unlock()
            }
        }
        else if( out != null ) {
            print out.toString()
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

    /**
     * Get the map of *input* objects by the given {@code InParam} type
     *
     * @param types One ore more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    def <T extends InParam> Map<T,Object> getInputsByType( Class<T>... types ) {

        def result = [:]
        inputs.findAll() { types.contains(it.key.class) }.each { result << it }
        return result
    }

    /**
     * Get the map of *output* objects by the given {@code InParam} type
     *
     * @param types One ore more subclass of {@code InParam}
     * @return An associative array containing all the objects for the specified type
     */
    def Map<OutParam,Object> getOutputsByType( Class... types ) {
        def result = [:]
        outputs.findAll() { types.contains(it.key.class) }.each { result << it }
        return result
    }

    /**
     * @return A map containing the task environment defined as input declaration by this task
     */
    Map<String,String> getInputEnvironment() {
        final Map<String,String> environment = [:]
        getInputsByType( EnvInParam ).each { param, value ->
            environment.put( param.name, value?.toString() )
        }
        return environment
    }

    /**
     * @return The location of the environment script used by this task
     */
    Path getCmdEnvironmentFile() { workDirectory.resolve('.command.env') }

    /**
     * @return The location of the task script script to be executed
     */
    Path getCmdScriptFile() { workDirectory.resolve('.command.sh') }

    /**
     * @return The location of the data file to be provided to this task
     */
    Path getCmdInputFile() { workDirectory.resolve('.command.in') }

    /**
     * @return The location of the data output by this task
     */
    Path getCmdOutputFile() { workDirectory.resolve('.command.out') }

    /**
     * @return The location of the file where the task exit status is saved
     */
    Path getCmdExitFile() { workDirectory.resolve('.exitcode') }

    /**
     * @return The location of the created when the task starts
     */
    Path getCmdStartedFile() { workDirectory.resolve('.command.started') }

    /**
     * @return The location of the wrapper BASH script user to launch the user target script
     */
    Path getCmdWrapperFile() { workDirectory.resolve('.command.run')  }

}

