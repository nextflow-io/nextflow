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
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Session

/**
 * The 'contract' of an abstract workflow command executor.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface TaskProcessor {

    /**
     * Define the execution environment of the command to be executed
     *
     * @param environment
     */
    TaskProcessor environment(Map<String,String> environment)

    /**
     * Defines the inputs dataflow
     *
     * @param inputs
     */
    TaskProcessor input(Map<String,?> inputs)

    /**
     * Defines the outputs dataflow
     */
    TaskProcessor output(Map<String,DataflowWriteChannel> outputs)

    TaskProcessor output(String... files)

    TaskProcessor echo( boolean value )

    TaskProcessor shareWorkDir( boolean value )

    TaskProcessor shell( String value )

    TaskProcessor validExitCodes( List<Integer> values )

    TaskProcessor errorStrategy( ErrorStrategy value )


    /**
     * Defines the maximum number of threads allowed to run this job concurrently
     *
     * @param max
     * @return
     */
    TaskProcessor threads( int max );

    /**
     * Define the processor 'name' attribute
     *
     * @param name
     * @return
     */
    TaskProcessor name( String name )

    TaskProcessor script( Closure closure )

    TaskProcessor script( String shell, Closure closure )

    /**
     * The code to be executed
     *
     * @param code
     */
    def run()

    /**
     * @return The current processing session
     */
    Session getSession()

    /**
     * Getter method for processor 'name' attribute
     * @return The processor 'name'
     */
    String getName()

    /**
     * Getter method for input channels
     * @param name The name of a input channel
     * @return The {@code DataflowReadChannel} instance of the requested channel or {@code null} if does not exist
     */
    DataflowReadChannel getInput( String name )

    /**
     * Getter method for output channel
     * @param name The name of a output channel
     * @return The {@code DataflowWriteChannel} instance of the requested channel or {@code null} if does not exist
     */
    DataflowWriteChannel getOutput( String name )

    /**
     * @return Whenever the task stdout console redirection is active
     */
    boolean getEcho()

    /**
     * @return The current environment map
     */
    Map<String,String > getEnvironment()

    /**
     * @return The maximum number of thread that can be used by the processor
     */
    int getThreads()

    /**
     * @return Whenever the same working direction have to be used for all tasks execution by this processor
     */
    boolean getShareWorkDir()

    String getShell()

    List<Integer> getValidExitCodes()

    ErrorStrategy getErrorStrategy()
}
