/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

package nextflow.trace

import es.bsc.cepbatools.extrae.Wrapper
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
/**
 * Profile observer producing an Extra tracing file
 * <p>
 * This class depends on Extrae 2.5.0 native library that need
 * to be compiled and installed in the hosting system.
 * <p>
 * When the pipeline starts the class looks for the Extrae
 * configuration file define by the env variable {@code EXTRAE_CONFIG_FILE}
 * if the config file is not available it will use the default
 * configuration file shipped with Nextflow.
 *
 *
 * @link http://www.bsc.es/computer-sciences/extrae
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Juan Gonzalez <juan.gonzalez@bsc.es>
 */
@Slf4j
@CompileStatic
class ExtraeTraceObserver implements TraceObserver {

    private Session session

    private Map<Integer,String> tasksMap = [:]

    @Override
    boolean enableMetrics() {
        return true
    }

    /**
     * Initialize the Extrae sub-system. It look the Extrae configuration
     * file defined by the variable {@code EXTRAE_CONFIG_FILE}.
     *
     * @throws IllegalArgumentException
     *      If the {@code EXTRAE_CONFIG_FILE} environment variable is not defined
     * @thows IOException
     *      If the configuration file cannot be written in the path specified by {@code EXTRAE_CONFIG_FILE}
     */
    @Override
    void onFlowStart(Session session) {

        String configVar = System.getenv('EXTRAE_CONFIG_FILE')
        log.debug "Starting tracing with Extrae - EXTRAE_CONFIG_FILE=$configVar"
        this.session = session

        if( !configVar )
            throw new IllegalArgumentException("Missing variable EXTRAE_CONFIG_FILE")

        def configFile = new File(configVar)
        if( !configFile.exists() )
            createConfigFile(configFile)

        Wrapper.defineEventType(1001, "Nextflow task submission", 0, null, null)

        Wrapper.resumeVirtualThread(1)
        Wrapper.Event(1000,1)
    }

    /**
     * This method is invoked when the flow is going to complete
     */
    @Override
    void onFlowComplete() {
        Wrapper.resumeVirtualThread(1)
        Wrapper.Event(1000,0)

        int size = tasksMap.size()
        def tasksTypes = new long[size]
        def tasksNames = new String [size]

        tasksMap.eachWithIndex{ int id, String name, int p ->
            tasksTypes[p] = id
            tasksNames[p] = name
        }

        Wrapper.defineEventType(1002, "Nextflow Process Types", size, tasksTypes, tasksNames)

        // reminder message for the user
        if( !session.isAborted() )
            log.info 'Hint: now you can generate the Paraver trace file by using the following command\n  ${EXTRAE_HOME}/bin/mpi2prv -task-view -f TRACE.mpits -o <file name.prv>\n'
    }

    /*
     * Invoked when the process is created.
     */
    void onProcessCreate( TaskProcessor process ) {
        // keep track of the process ID and name
        tasksMap.put( process.id, process.name )
    }


    @Override
    void onProcessSubmit(TaskHandler handler) {

    }

    @Override
    void onProcessStart( TaskHandler handler ) {

        long taskId = handler.getTask().id as long
        long procId = handler.getTask().processor.id

        Wrapper.resumeVirtualThread( taskId+1 )
        Wrapper.Event(1001, taskId )    // Task ID
        Wrapper.Event(1002, procId )    // Process ID
        Wrapper.Event(1003, 0)          // % cpu
        Wrapper.Event(1004, 0)          // rss
        Wrapper.Event(1005, 0)          // virtual memory
        Wrapper.Event(1006, 0)          // read chars
        Wrapper.Event(1007, 0)          // write chars
        Wrapper.Event(1008, 0)          // number of read ops
        Wrapper.Event(1009, 0)          // number of write ops

        Wrapper.suspendVirtualThread()

    }

    @Override
    void onProcessComplete( TaskHandler handler ) {

        final taskId = handler.getTask().id as long
        final record = handler.getTraceRecord()

        Wrapper.resumeVirtualThread( taskId+1 )
        Wrapper.Event(1001, 0)
        Wrapper.Event(1002, 0)
        Wrapper.Event(1003, (record.get('%cpu') as Float ?: 0f) * 10L as long)      // note: since the cpu usage percentage can contain a decimal digit, the value is multiply by 10
        Wrapper.Event(1004, (record.get('rss') ?: 0) as long)                       // rss (bytes)
        Wrapper.Event(1005, (record.get('vmem') ?: 0) as long)                      // virtual memory (bytes)
        Wrapper.Event(1006, (record.get('rchar') ?: 0) as long)                     // read chars (bytes)
        Wrapper.Event(1007, (record.get('wchar') ?: 0) as long)                     // write chars (bytes)
        Wrapper.Event(1008, (record.get('syscr') ?: 0) as long)                     // number of read ops
        Wrapper.Event(1009, (record.get('syscw') ?: 0) as long)                     // number of write ops

        Wrapper.suspendVirtualThread()

    }

    void onProcessCached( TaskHandler handler ) {
        //onProcessComplete(handler)
    }

    /**
     * Extract the configuration file embedded in the distribution JAR and copy into the
     * specified file path
     *
     * @param configFile The file path where the config file has to be created
     */
    @PackageScope
    void createConfigFile(File configFile) {

        if( configFile.parent && !configFile.parentFile.exists() && !configFile.parentFile.mkdirs() )
            throw new IOException("Cannot create Extrae config file: $configFile -- Check write file system permissions")

        this.getClass()
                .getResourceAsStream("/extrae.xml")
                .withReader { source ->
                    configFile.withOutputStream { target ->
                        target << source
                    }
        }
    }
}
