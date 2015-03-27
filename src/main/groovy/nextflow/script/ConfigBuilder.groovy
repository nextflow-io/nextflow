/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
import static nextflow.Const.APP_HOME_DIR
import static nextflow.util.ConfigHelper.parseValue

import java.nio.file.Path
import java.nio.file.Paths

import groovy.util.logging.Slf4j
import nextflow.cli.CliOptions
import nextflow.cli.CmdNode
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.util.HistoryFile

/**
 * Builds up the Nextflow configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigBuilder {

    CliOptions options

    CmdRun cmdRun

    CmdNode cmdNode

    Path baseDir

    Path workDir = Paths.get('.').complete()

    ConfigBuilder setOptions( CliOptions options ) {
        this.options = options
        return this
    }

    ConfigBuilder setCmdRun( CmdRun cmdRun ) {
        this.cmdRun = cmdRun
        return this
    }

    ConfigBuilder setBaseDir( Path file ) {
        this.baseDir = file
        return this
    }

    ConfigBuilder setCmdNode( CmdNode node ) {
        this.cmdNode = node
        return this
    }


    static private wrapValue( value ) {
        if( !value )
            return ''

        value = value.toString().trim()
        if( value == 'true' || value == 'false')
            return value

        if( value.isNumber() )
            return value

        return "'$value'"
    }

    /**
     * Transform the specified list of string to a list of files, verifying their existence.
     * <p>
     *     If a file in the list does not exist an exception of type {@code CliArgumentException} is thrown.
     * <p>
     *     If the specified list is empty it tries to return of default configuration files located at:
     *     <li>$HOME/.nextflow/taskConfig
     *     <li>$PWD/nextflow.taskConfig
     *
     * @param files
     * @return
     */
    def List<Path> validateConfigFiles( List<String> files ) {

        def result = []
        if ( files ) {
            files.each { String fileName ->
                def thisFile = Paths.get(fileName)
                if(!thisFile.exists()) {
                    throw new AbortOperationException("The specified configuration file does not exist: $thisFile -- check the name or choose another file")
                }
                result << thisFile
            }
            return result
        }

        /*
         * config file in the nextflow home
         */
        def home = APP_HOME_DIR.resolve('config')
        if( home.exists() ) {
            log.debug "Found config home: $home"
            result << home
        }

        /**
         * Config file in the pipeline base dir
         */
        def base = null
        if( baseDir && baseDir.complete() != Paths.get('.').complete() ) {
            base = baseDir.resolve('nextflow.config')
            if( base.exists() ) {
                log.debug "Found config base: $base"
                result << base
            }
        }

        /**
         * Local or user provided file
         */
        def local = Paths.get('nextflow.config')
        if( local.exists() && local != base ) {
            log.debug "Found config local: $local"
            result << local
        }

        def customConfigs = []
        if( options?.userConfig ) customConfigs.addAll(options.userConfig)
        if( cmdRun?.runConfig ) customConfigs.addAll(cmdRun.runConfig)
        if( customConfigs ) {
            for( String item : customConfigs ) {
                def configFile = Paths.get(item)
                if(!configFile.exists()) {
                    throw new AbortOperationException("The specified configuration file does not exist: $configFile -- check the name or choose another file")
                }

                log.debug "User config file: $configFile"
                result << configFile
            }
        }

        return result
    }

    def Map buildConfig( List<Path> files ) {

        final Map<String,String> vars = cmdRun?.env
        final boolean exportSysEnv = cmdRun?.exportSysEnv

        def items = []
        files?.each { Path file ->
            log.debug "Parsing config file: ${file.complete()}"
            if (!file.exists()) {
                log.warn "The specified configuration file cannot be found: $file"
            }
            else {
                items << file
            }
        }

        Map env = [:]
        if( exportSysEnv ) {
            log.debug "Adding current system environment to session environment"
            env.putAll(System.getenv())
        }
        if( vars ) {
            log.debug "Adding following variables to session environment: $vars"
            env.putAll(vars)
        }

        // set the cluster options for the node command
        if( cmdNode?.clusterOptions )  {
            def str = new StringBuilder()
            cmdNode.clusterOptions.each { k, v ->
                str << "cluster." << k << '=' << wrapValue(v) << '\n'
            }
            items << str
        }

        // -- add the executor obj from the command line args
        if( cmdRun?.clusterOptions )  {
            def str = new StringBuilder()
            cmdRun.clusterOptions.each { k, v ->
                str << "cluster." << k << '=' << wrapValue(v) << '\n'
            }
            items << str
        }

        if( cmdRun?.executorOptions )  {
            def str = new StringBuilder()
            cmdRun.executorOptions.each { k, v ->
                str << "executor." << k << '=' << wrapValue(v) << '\n'
            }
            items << str
        }

        buildConfig0( env, items )
    }


    def Map buildConfig0( Map env, List configEntries )  {
        assert env != null

        ConfigObject result = new ConfigSlurper().parse('env{}; session{}; params{}; process{}; executor{} ')

        // add the user specified environment to the session env
        env.sort().each { name, value -> result.env.put(name,value) }

        if( configEntries ) {
            // the configuration object binds always the current environment
            // so that in the configuration file may be referenced any variable
            // in the current environment
            final binding = new HashMap(System.getenv())
            binding.putAll(env)
            binding.put('workDir', workDir)
            binding.put('baseDir', baseDir)

            configEntries.each { entry ->

                def text = entry instanceof Path ? entry.text : entry.toString()

                if ( text ) {
                    def cfg = new ConfigSlurper()
                    cfg.setBinding(binding)
                    try {
                        result.merge( cfg.parse(text) )
                    }
                    catch( Exception e ) {
                        def message = (entry instanceof Path ? "Unable to parse config file: '$entry' " : "Unable to parse configuration ")
                        throw new ConfigParseException(message,e)
                    }
                }
            }

        }

        return result
    }

    private String normalizeResumeId( String uniqueId ) {
        if( !uniqueId )
            return null

        if( uniqueId == 'last' ) {
            uniqueId = HistoryFile.history.retrieveLastUniqueId()
            if( !uniqueId ) {
                log.warn "It appears you have never run this pipeline before -- Option `-resume` is ignored"
            }
            return uniqueId
        }

        if( HistoryFile.history.findUniqueId(uniqueId))
            return uniqueId

        throw new AbortOperationException("Can't find a run with the specified id: ${uniqueId} -- Execution can't be resumed")
    }

    private configRunOptions(Map config) {

        // -- override 'process' parameters defined on the cmd line
        cmdRun.process.each { name, value ->
            config.process[name] = parseValue(value)
        }

        // -- check for the 'continue' flag
        if( cmdRun.resume ) {
            config.session.uniqueId = normalizeResumeId(cmdRun.resume)
        }
        // -- other configuration parameters
        if( cmdRun.poolSize ) {
            config.poolSize = cmdRun.poolSize
        }
        if( cmdRun.queueSize ) {
            config.executor.queueSize = cmdRun.queueSize
        }
        if( cmdRun.pollInterval ) {
            config.executor.pollInterval = cmdRun.pollInterval
        }

        // -- add the command line parameters to the 'taskConfig' object
        cmdRun.params?.each { name, value ->
            config.params.put(name, parseValue(value))
        }

        if( cmdRun.withoutDocker && config.docker instanceof Map ) {
            // disable docker execution
            log.debug "Disabling execution in Docker contained as request by cli option `-without-docker`"
            config.docker.enabled = false
        }

        if( cmdRun.withDocker ) {
            log.debug "Enabling execution in Docker container as request by cli option `-with-docker ${cmdRun.withDocker}`"
            if( !config.containsKey('docker') )
                config.docker = [:]
            if( !(config.docker instanceof Map) )
                throw new AbortOperationException("Invalid `docker` definition in the config file")

            config.docker.enabled = true
            if( cmdRun.withDocker != '-' ) {
                // this is supposed to be a docker image name
                config.process.container = cmdRun.withDocker
            }
            else if( config.docker.image ) {
                config.process.container = config.docker.image
            }

            if( !hasContainerDirective(config.process) )
                throw new AbortOperationException("You request to run with Docker but no image has been specified")

        }
    }

    /**
     * Verify that configuration for process contains at last one `container` directive
     *
     * @param process
     * @return {@code true} when a `container` is defined or {@code false} otherwise
     */
    protected boolean hasContainerDirective(process)  {

        if( process instanceof Map ) {
            if( process.container )
                return true

            def result = process
                            .findAll { String name, value -> name.startsWith('$') && value instanceof Map }
                            .find { String name, Map value -> value.container as boolean }  // the first non-empty `container` string

            return result as boolean
        }

        return false
    }



    /**
     * Given the command line options {@code CliOptions} object
     * read the application configuration file and returns the
     * config object
     *
     * @param options The {@code CliOptions} as specified by the user
     * @return A the application options hold in a {@code Map} object
     */
    ConfigObject build() {

        // -- configuration file(s)
        def configFiles = validateConfigFiles(options?.config)
        def config = buildConfig(configFiles)

        if( cmdRun )
            configRunOptions(config)

        // convert the ConfigObject to plain map
        // this because when accessing a non-existing entry in a ConfigObject it return and empty map as value
        return config
    }

}
