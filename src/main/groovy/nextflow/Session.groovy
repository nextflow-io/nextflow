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

package nextflow
import com.google.common.collect.LinkedHashMultimap
import com.google.common.collect.Multimap
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.processor.LocalScriptProcessor
import nextflow.processor.NopeScriptProcessor
import nextflow.processor.OgeScriptProcessor
import nextflow.processor.Processor
import nextflow.processor.TaskDef
import nextflow.script.AbstractScript

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Session {


    /**
     * The class to be used create a {@code Processor} instance
     */
    Class<? extends Processor> processorClass = LocalScriptProcessor.class

    /**
     * Keep a list of all processor created
     */
    List<DataflowProcessor> allProcessors = []

    /**
     * Keep the list all executed tasks
     * note: LinkedHashMultimap preserves insertion order of entries, as well as the insertion order of keys, and the set of values associated with any one key.
     */
    Multimap<Processor, TaskDef> tasks = LinkedHashMultimap.create()

    /**
     * The directory where tasks create their execution folders and store results, bye default the current directory
     */
    File workDirectory = new File('.')


    /**
     * Holds the configuration object
     */
    def Map config

    /**
     * Creates a new session with an 'empty' (default) configuration
     */
    def Session() {
        this([:])
    }


    /**
     * Creates a new session using the configuration properties provided
     *
     * @param config
     */
    def Session( Map config ) {
        assert config != null
        this.config = config

        // normalize config object
        if( config.task == null ) config.task = [:]
        if( config.env == null ) config.env = [:]

        this.processorClass = loadProcessorClass(config.task.processor?.toString())
    }


    protected Class<? extends Processor> loadProcessorClass(String processorType) {

        def className
        if ( !processorType ) {
            className = LocalScriptProcessor.name
        }
        else if ( processorType.toLowerCase() == 'local' ) {
            className = LocalScriptProcessor.name
        }
        else if ( processorType.toLowerCase() in ['sge','oge'] ) {
            className = OgeScriptProcessor.name
        }
        else if ( processorType.toLowerCase() == 'nope' ) {
            className = NopeScriptProcessor.name
        }
        else {
            className = processorType
        }

        log.debug "Loading processor class: ${className}"
        try {
            Thread.currentThread().getContextClassLoader().loadClass(className) as Class<Processor>
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Cannot find a valid class for specified processor type: '${processorType}'")
        }

    }


    /**
     * Create an instance of the task {@code Processor}
     * @return
     */
    Processor createProcessor(AbstractScript script = null, boolean bindOnTermination = false) {

        // -- create a new processor instance
        def processor = processorClass.newInstance( this, script, bindOnTermination )

        // -- inject attributes defined by the 'config.task' element
        if ( config.task instanceof Map ) {

            def methods = processor.metaClass.getMethods().findAll{ MetaMethod m -> m.isPublic() && m.getParameterTypes().size()==1 }
            def names = methods *. getName()

            config.task.each { String key, Object value ->

                def i = names.indexOf(key)
                if ( i != -1 ) {

                    try {
                        methods[i].invoke(processor,value)
                    }
                    catch( Exception e ) {
                        def mType = methods[i].getNativeParameterTypes()[0]?.simpleName
                        def vType = value?.class?.simpleName
                        log.warn "Task attribute '$key' requires a value of type: '${mType}' -- entered value: '${value} of type: '${vType}'" , e
                    }
                }
            }
        }


        return processor

    }


    /**
     * Await the termination of all processors
     */
    void awaitTermination() {
        allProcessors *. join()
    }

//    /**
//     * Create a table report of all executed or running tasks
//     *
//     * @return A string table formatted displaying the tasks information
//     */
//    String tasksReport() {
//
//        TableBuilder table = new TableBuilder()
//                .head('name')
//                .head('id')
//                .head('status')
//                .head('path')
//                .head('exit')
//
//        tasks.entries().each { Map.Entry<Processor, TaskDef> entry ->
//            table << entry.key.name
//            table << entry.value.id
//            table << entry.value.status
//            table << entry.value.workDirectory
//            table << entry.value.exitCode
//            table << table.closeRow()
//        }
//
//        table.toString()
//
//    }

}
