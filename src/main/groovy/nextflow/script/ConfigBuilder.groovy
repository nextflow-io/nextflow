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

import static nextflow.util.ConfigHelper.parseValue

import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.ExitCode
import nextflow.cli.CliOptions
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.util.HistoryFile

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigBuilder {

    CliOptions options

    CmdRun runOptions


    ConfigBuilder setOptions( CliOptions options ) {
        this.options = options
        return this
    }

    ConfigBuilder setCmdRun( CmdRun cmdRun ) {
        this.runOptions = cmdRun
        return this
    }


    /**
     * Converts a {@code ConfigObject} to a plain {@code Map}
     *
     * @param config
     * @return
     */
    static Map configToMap( ConfigObject config ) {
        assert config != null

        Map result = new LinkedHashMap(config.size())
        config.keySet().each { name ->
            def value = config.get(name)
            if( value instanceof ConfigObject ) {
                result.put( name, configToMap(value))
            }
            else if ( value != null ){
                result.put( name, value )
            }
            else {
                result.put( name, null )
            }
        }

        return result
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
    def List<File> validateConfigFiles( List<String> files ) {

        def result = []
        if ( files ) {
            files.each { String fileName ->
                def thisFile = new File(fileName)
                if(!thisFile.exists()) {
                    throw new AbortOperationException("The specified configuration file does not exist: $thisFile -- check the name or choose another file")
                }
                result << thisFile
            }
            return result
        }

        def home = new File(Const.APP_HOME_DIR, 'config')
        if( home.exists() ) result << home

        def local = new File('nextflow.config')
        if( local.exists() ) result << local

        return result
    }


    def Map buildConfig( List<File> files ) {

        final Map<String,String> vars = runOptions?.env
        final boolean exportSysEnv = runOptions?.exportSysEnv

        def texts = []
        files?.each { File file ->
            log.debug "Parsing config file: ${file.absoluteFile}"
            if (!file.exists()) {
                log.warn "The specified configuration file cannot be found: $file"
            }
            else {
                texts << file.text
            }
        }

        Map<String,String> env = [:]
        if( exportSysEnv ) {
            log.debug "Adding current system environment to session environment"
            env.putAll(System.getenv())
        }
        if( vars ) {
            log.debug "Adding following variables to session environment: $vars"
            env.putAll(vars)
        }

        // -- add the daemon obj from the command line args
        if( options.daemonOptions )  {
            def str = new StringBuilder()
            options.daemonOptions.each { k, v ->
                str << "daemon." << k << '=' << wrapValue(v) << '\n'
            }
            texts << str.toString()
        }

        if( runOptions?.executorOptions )  {
            def str = new StringBuilder()
            runOptions.executorOptions.each { k, v ->
                str << "executor." << k << '=' << wrapValue(v) << '\n'
            }
            texts << str.toString()
        }

        buildConfig0( env, texts )
    }


    def static Map buildConfig0( Map env, List<String> confText )  {
        assert env != null

        ConfigObject result = new ConfigSlurper().parse('env{}; session{}; params{}; process{}; executor{} ')

        // add the user specified environment to the session env
        env.sort().each { name, value -> result.env.put(name,value) }

        if( confText ) {
            // the configuration object binds always the current environment
            // so that in the configuration file may be referenced any variable
            // in the current environment
            final binding = new HashMap(System.getenv())
            binding.putAll(env)

            confText.each { String text ->
                if ( text ) {
                    def cfg = new ConfigSlurper()
                    cfg.setBinding(binding)
                    try {
                        result.merge( cfg.parse(text) )
                    }
                    catch( Exception e ) {
                        throw new ConfigParseException("Failed to parse config file",e)
                    }
                }
            }

        }

        // convert the ConfigObject to plain map
        // this because when accessing a non-existing entry in a ConfigObject it return and empty map as value
        return configToMap( result )
    }

    private configRunOptions(Map config) {

        // -- override 'process' parameters defined on the cmd line
        runOptions.process.each { name, value ->
            config.process[name] = parseValue(value)
        }

        // -- check for the 'continue' flag
        if( runOptions.resume ) {
            def uniqueId = runOptions.resume
            if( uniqueId == 'last' ) {
                uniqueId = HistoryFile.history.retrieveLastUniqueId()
                if( !uniqueId ) {
                    log.error "It appears you have never executed it before -- Cannot use the '-resume' command line option"
                    System.exit(ExitCode.MISSING_UNIQUE_ID)
                }
            }
            config.session.uniqueId = uniqueId
        }

        // -- other configuration parameters
        if( runOptions.poolSize ) {
            config.poolSize = runOptions.poolSize
        }
        if( runOptions.queueSize ) {
            config.executor.queueSize = runOptions.queueSize
        }
        if( runOptions.pollInterval ) {
            config.executor.pollInterval = runOptions.pollInterval
        }

        // -- add the command line parameters to the 'taskConfig' object
        runOptions.params?.each { name, value ->
            config.params.put(name, parseValue(value))
        }
    }

    /**
     * Given the command line options {@code CliOptions} object
     * read the application configuration file and returns the
     * config object
     *
     * @param options The {@code CliOptions} as specified by the user
     * @return A the application options hold in a {@code Map} object
     */
    Map build() {

        // -- configuration file(s)
        def configFiles = validateConfigFiles(options?.config)
        def config = buildConfig(configFiles)

        if( runOptions )
            configRunOptions(config)

        return config
    }

}
