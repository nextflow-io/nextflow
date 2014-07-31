/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.ExitCode
import nextflow.daemon.DaemonLauncher
import nextflow.exception.MissingScriptException
import nextflow.exception.ConfigParseException
import nextflow.executor.ServiceName
import nextflow.util.FileHelper
import nextflow.util.LoggerHelper
import nextflow.util.ServiceDiscover

/**
 * Main application entry point. It parses the command line and
 * launch the pipeline execution.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Launcher {

    /**
     * Create the application command line parser
     *
     * @return An instance of {@code CliBuilder}
     */

    static JCommander jcommander

    static boolean fullVersion

    static CliOptions options

    static CmdX command

    @PackageScope
    static void parseMainArgs(String... args) {

        def commands = [
                clone:   new CmdClone(),
                history: new CmdHistory(),
                info:    new CmdInfo(),
                list:    new CmdList(),
                pull:    new CmdPull(),
                run:     new CmdRun(),
                drop:    new CmdDrop(),
                search:  new CmdSearch()
        ]

        def normalizedArgs = CliOptions.normalizeArgs(args)
        options = new CliOptions()
        jcommander = new JCommander(options)
        commands.each { name, cmd -> jcommander.addCommand(name, cmd) }
        jcommander.setProgramName( Const.APP_NAME )
        jcommander.parse( normalizedArgs as String[] )
        fullVersion = '-version' in normalizedArgs
        command = commands.get( jcommander.getParsedCommand() )
        if( command )
            command.setOptions(options)
    }


    /**
     * Program entry method
     *
     * @param args The program options as specified by the user on the CLI
     */
    public static void main(String... args)  {

        def scriptBaseName = null
        try {
            // -- parse the program arguments - and - configure the logger
            parseMainArgs(args)
            LoggerHelper.configureLogger(options)
            log.debug "${Const.APP_NAME} ${args?.join(' ')}"

            // -- print out the version number, then exit
            if ( options.version ) {
                println getVersion(fullVersion)
                System.exit(ExitCode.OK)
            }

            // -- print out the program help, then exit
            if( options.help || (!command && !options.isDaemon())) {
                println Const.LOGO
                jcommander.usage()
                System.exit(ExitCode.OK)
            }

            if( !options.quiet && !(command instanceof CmdInfo) ) {
                println "N E X T F L O W  ~  version ${Const.APP_VER}"
            }

            // -- launch daemon
            if( options.isDaemon() ) {
                log.debug "Launching cluster daemon"
                launchDaemon()
                return
            }

            // launch the command
            command.run()

        }

        catch( ParameterException e ) {
            // print command line parsing errors
            // note: use system.err.println since if an exception is raised
            //       parsing the cli params the logging is not configured
            System.err.println "${e.getMessage()} -- Check the available command line parameters and syntax using '-h'"
            System.exit( ExitCode.INVALID_COMMAND_LINE_PARAMETER )
        }

        catch (MissingScriptException e ) {
            log.error "The specified script does not exist: '${e.scriptFile}'\n", e
            System.exit( ExitCode.MISSING_SCRIPT_FILE )
        }

        catch( ConfigParseException e )  {
            log.error "${e.message}\n${e.cause}\n\n"
            System.exit(ExitCode.INVALID_CONFIG)
        }

        catch( Throwable fail ) {
            log.error("${fail.toString()} -- See the file '.nextflow.log' for more error details", fail)
            System.exit( ExitCode.UNKNOWN_ERROR )
        }

    }


    /**
     * Print the application version number
     * @param full When {@code true} prints full version number including build timestamp
     * @return The version number string
     */
    static String getVersion(boolean full = false) {

        if ( full ) {
            Const.LOGO
        }
        else {
            "${Const.getAPP_NAME()} version ${Const.APP_VER}.${Const.APP_BUILDNUM}"
        }

    }

    /**
     * Launch the daemon service
     *
     * @param config The nextflow configuration map
     */
    static launchDaemon() {

        // -- check file system providers
        FileHelper.checkFileSystemProviders()

        // create the config object
        def config = new ConfigBuilder().setOptions(options).build()

        def daemonConfig = config.daemon instanceof Map ? config.daemon : [:]
        log.debug "Daemon config > $daemonConfig"


        DaemonLauncher instance
        def name = daemonConfig.name as String
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
        instance.launch(daemonConfig)
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
        for( Class item : ServiceDiscover.load(DaemonLauncher).iterator() ) {
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
            throw new IllegalStateException("No daemon services are available -- Cannot launch Nextflow in damon mode")

        return loader.next()
    }

}
