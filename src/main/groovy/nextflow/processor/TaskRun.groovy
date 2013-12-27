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
import java.nio.file.Paths
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.script.ScriptType

/**
 * Models a task instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ToString( includePackage = false, includeNames = true, includes = 'id,index,name,status,exitCode' )
class TaskRun {

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
    int exitStatus = Integer.MAX_VALUE

    /**
     * Flag set when the bind stage is completed successfully
     */
    boolean canBind

    /**
     * The return std out
     */
    def stdout

    /**
     * @return The task produced stdout result as string
     */
    def String getStdout() {

        if( stdout instanceof Path ) {
            return stdout.exists() ? stdout.text : null
        }
        else if( stdout != null ) {
            return stdout.toString()
        }

        return null
    }

    /**
     * Print to the current system console the task produced stdout
     */
    void echoStdout() {

        // print the stdout
        if( stdout instanceof Path ) {
            if( !stdout.exists() ) {
                log.debug "Echo file does not exist: ${out}"
                return
            }

            echoLock.lock()
            try {
                stdout.withReader {  System.out << it }
            }
            finally {
                echoLock.unlock()
            }
        }
        else if( stdout != null ) {
            print stdout.toString()
        }

    }


    /**
     * The directory used to run the task
     */
    Path workDirectory

    /**
     * The type of the task code: native or external system command
     */
    ScriptType type

    /**
     * The runtime exception when execution groovy native code
     */
    Throwable error

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


    static Path CMD_ENV = Paths.get('.command.env')
    static Path CMD_SCRIPT = Paths.get('.command.sh')
    static Path CMD_INFILE = Paths.get('.command.in')
    static Path CMD_OUTFILE = Paths.get('.command.out')
    static Path CMD_EXIT = Paths.get('.exitcode')
    static Path CMD_START = Paths.get('.command.started')
    static Path CMD_RUN = Paths.get('.command.run')
    static Path CMD_CONTEXT = Paths.get('.command.ctx')

    /**
     * @return The location of the environment script used by this task
     */
    Path getCmdEnvironmentFile() { workDirectory.resolve(CMD_ENV) }

    /**
     * @return The location of the task script script to be executed
     */
    Path getCmdScriptFile() { workDirectory.resolve(CMD_SCRIPT) }

    /**
     * @return The location of the data file to be provided to this task
     */
    Path getCmdInputFile() { workDirectory.resolve(CMD_INFILE) }

    /**
     * @return The location of the data output by this task
     */
    Path getCmdOutputFile() { workDirectory.resolve(CMD_OUTFILE) }

    /**
     * @return The location of the file where the task exit status is saved
     */
    Path getCmdExitFile() { workDirectory.resolve(CMD_EXIT) }

    /**
     * @return The location of the created when the task starts
     */
    Path getCmdStartedFile() { workDirectory.resolve(CMD_START) }

    /**
     * @return The location of the wrapper BASH script user to launch the user target script
     */
    Path getCmdWrapperFile() { workDirectory.resolve(CMD_RUN)  }

    /**
     * @return The location of the file that holds the cached process context map
     */
    Path getCmdContextFile() { workDirectory.resolve(CMD_CONTEXT) }

}

