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

package nextflow.script
import groovy.transform.PackageScope
import nextflow.Global
import nextflow.Session
import nextflow.processor.ProcessFactory
import nextflow.processor.TaskProcessor
/**
 * Any user defined script will extends this class, it provides the base execution context
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class BaseScript extends Script {

    protected BaseScript() { }

    protected BaseScript(Binding binding) {
        super(binding)
    }

    /**
     * This method is get invoked by the DSL parser
     * @param processNames
     */
    protected void init( List<String> processNames ) {
        this.processNames = processNames
    }

    /*
     * The script execution session, declare it private to prevent the user script to be able to access it
     */
    private Session sessionObj

    /**
     * The list of process defined in the pipeline script
     */
    private List<String> processNames

    @PackageScope
    void setSession( Session value )  {
        sessionObj = value
    }

    private Session getSession() {
        if( !sessionObj ) {
            sessionObj = Global.session as Session
        }
        return sessionObj
    }

    @PackageScope
    List<String> getProcessNames() {
        processNames
    }

    /**
     * Holds the configuration object which will used to execution the user tasks
     */
    Map getConfig() { getSession().getConfig() }

    @Lazy
    InputStream stdin = { System.in }()

    private TaskProcessor taskProcessor

    /** Access to the last *process* object -- only for testing purpose */
    @PackageScope
    TaskProcessor getTaskProcessor() { taskProcessor }

    private result

    /** Access to the last *process* result -- only for testing purpose */
    @PackageScope
    Object getResult() { result }

    private ProcessFactory processFactory

    private getProcessFactory() {
        if( !processFactory ) {
            processFactory = new ProcessFactory(this, getSession())
        }
        return processFactory
    }

    @PackageScope
    void setProcessFactory( ProcessFactory value ) {
        this.processFactory = value
    }

    /**
     * Enable disable task 'echo' configuration property
     * @param value
     */
    def void echo(boolean value = true) {
        config.process.echo = value
    }


    /**
     * Method to which is mapped the *process* declaration when using the following syntax:
     *    <pre>
     *        process myProcess( option: value[, [..]] ) &#123;
     *        ..
     *        &#125;
     *    </pre>
     *
     * @param args The process options specified in the parenthesis after the process name, as shown in the above example
     * @param name The name of the process, e.g. {@code myProcess} in the above example
     * @param body The body of the process declaration. It holds all the process definitions: inputs, outputs, code, etc.
     */
    protected process( Map<String,?> args, String name, Closure body ) {
        // create the process
        taskProcessor = getProcessFactory().createProcessor(name, body, args)
        // launch it
        result = taskProcessor.run()
    }

    /**
     * Method to which is mapped the *process* declaration.
     *
     *    <pre>
     *        process myProcess( option: value[, [..]] ) &#123;
     *        ..
     *        &#125;
     *    </pre>
     * @param name The name of the process, e.g. {@code myProcess} in the above example
     * @param body The body of the process declaration. It holds all the process definitions: inputs, outputs, code, etc.
     */
    protected process( String name, Closure body ) {
        // create the process
        taskProcessor = getProcessFactory().createProcessor(name, body)
        // launch it
        result = taskProcessor.run()
    }

}
