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
import groovy.text.GStringTemplateEngine
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.processor.LocalScriptProcessor
import nextflow.processor.Processor
import nextflow.processor.TaskDef
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
     * The default environment used to run the tasks process
     */
    @Lazy
    Map<String,String> env = { createEnvMap() } ()


    private Map<String, String> createEnvMap() {

        // by default load the current environment variables
        def result = new HashMap<String,String>(System.getenv())


        // override current environment by file in the 'nextflow' home directory
        def home = new File(Const.APP_HOME_DIR,'environment')
        if ( home.exists() ) {
            log.debug "Loading home environment file: $home"
            result.putAll( parseEnvironmentFile(home, result) )
        }

        // override current environment by the file in the local directory
        def local = new File('.environment')
        if ( local.exists() ) {
            log.debug "Loading local environment file: $home"
            result.putAll( parseEnvironmentFile(home, result) )
        }


        return result
    }


    @TupleConstructor
    static class SoftReplaceMap implements Map {

        @Delegate
        Map map

        @Override
        def get(Object name ) {
            if( map.containsKey(name) ) return map[name]
            // fallback on Java system properties
            else if ( System.properties.containsKey(name) ) return System.properties[name]
            // otherwise return '' (an empty string)
            else return ''
        }
    }

    private Map<String,String> parseEnvironmentFile( File template, Map<String,String> binding ) {
        assert template
        assert binding != null


        def engine = new GStringTemplateEngine()
        def result = new LinkedHashMap<String,String>(binding)
        int c = 0
        template.text.readLines().each { line->

            def matcher = ( line =~~ /^\s*(\S+)\s*=\s*(.*)\s*$/ )
            if ( !matcher.matches()  )  {
                return
            }

            String name = matcher[0][1]
            String value = matcher[0][2]

            // when the quote is single-quote delimited
            // it is not interpreted
            matcher = ( value =~~ /'(.*)'/ )
            if( matcher.matches() ) {
                value = matcher[0][1]
                result.put(name,value)
            }
            // otherwise replace any variables
            else {
                matcher = ( value =~~ /"(.*)"/ )
                if( matcher.matches() ) {
                    value = matcher[0][1]
                }

                value = engine.createTemplate(value?.toString()).make(new SoftReplaceMap(result))
                result.put( name, value?.toString() )
            }

        }

        return result
    }


    /**
     * Create an instance of the task {@code Processor}
     * @return
     */
    Processor createProcessor(boolean bindOnTermination = false) {
        processorClass.newInstance( this, bindOnTermination )
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
