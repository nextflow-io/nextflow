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

package nextflow.cli
import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.config.ConfigBuilder
import nextflow.daemon.DaemonLauncher
import nextflow.util.ServiceName
import nextflow.util.ServiceDiscover
import picocli.CommandLine

/**
 * CLI-command NODE
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters
@CommandLine.Command
class CmdNode extends CmdBase {

    static final public NAME = 'node'

    @Override
    final String getName() { NAME }

    @DynamicParameter(names ='-cluster.', description='Define cluster config options')
    Map<String,String> clusterOptions = [:]

    @Parameter(names = ['-bg'], arity = 0, description = 'Start the cluster node daemon in background')
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @Parameter
    List<String> provider

    @Override
    void run() {
        launchDaemon(provider ? provider[0] : null)
    }


    /**
     * Launch the daemon service
     *
     * @param config The nextflow configuration map
     */
    protected launchDaemon(String name = null) {

        // create the config object
        def config = new ConfigBuilder()
                            .setOptions(launcher.options)
                            .setCmdNode(this)
                            .build()

        DaemonLauncher instance
        if( name ) {
            if( name.contains('.') ) {
                instance = loadDaemonByClass(name)
            }
            else {
                instance = loadDaemonByName(name)
            }
        }
        else {
            instance = loadDaemonFirst()
        }


        // launch it
        instance.launch(config)
    }

    /**
     * Load a {@code DaemonLauncher} instance of the its *friendly* name i.e. the name provided
     * by using the {@code ServiceName} annotation on the daemon class definition
     *
     * @param name The executor name e.g. {@code gridgain}
     * @return The daemon launcher instance
     * @throws IllegalStateException if the class does not exist or it cannot be instantiated
     */
    static DaemonLauncher loadDaemonByName( String name ) {

        Class<DaemonLauncher> clazz = null
        final itr = ServiceDiscover.load(DaemonLauncher).iterator()
        while( itr.hasNext() ) {
            final item = itr.next()
            log.debug "Discovered daemon class: ${item.name}"
            ServiceName annotation = item.getAnnotation(ServiceName)
            if( annotation && annotation.value() == name ) {
                clazz = item
                break
            }
        }

        if( !clazz )
            throw new IllegalStateException("Unknown daemon name: $name")

        try {
            clazz.newInstance()
        }
        catch( Exception e ) {
            throw new IllegalStateException("Unable to launch executor: $name", e)
        }
    }

    /**
     * Load a class implementing the {@code DaemonLauncher} interface by the specified class name
     *
     * @param name The fully qualified class name e.g. {@code nextflow.executor.LocalExecutor}
     * @return The daemon launcher instance
     * @throws IllegalStateException if the class does not exist or it cannot be instantiated
     */
    static DaemonLauncher loadDaemonByClass( String name ) {
        try {
            return (DaemonLauncher)Class.forName(name).newInstance()
        }
        catch( Exception e ) {
            throw new IllegalStateException("Cannot load daemon: ${name}")
        }
    }

    /**
     * @return The first available instance of a class implementing {@code DaemonLauncher}
     * @throws IllegalStateException when no class implementing {@code DaemonLauncher} is available
     */
    static DaemonLauncher loadDaemonFirst() {
        def loader = ServiceLoader.load(DaemonLauncher).iterator()
        if( !loader.hasNext() )
            throw new IllegalStateException("No cluster services are available -- Cannot launch Nextflow in cluster mode")

        return loader.next()
    }

}
