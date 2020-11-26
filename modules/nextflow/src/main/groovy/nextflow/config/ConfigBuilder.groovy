/*
 * Copyright 2020, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.config

import static nextflow.util.ConfigHelper.*

import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.cli.CliOptions
import nextflow.cli.CmdConfig
import nextflow.cli.CmdNode
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.trace.GraphObserver
import nextflow.trace.ReportObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
import nextflow.trace.WebLogObserver
import nextflow.util.HistoryFile
/**
 * Builds up the Nextflow configuration object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ConfigBuilder {

    static final public String DEFAULT_PROFILE = 'standard'

    CliOptions options

    CmdRun cmdRun

    CmdNode cmdNode

    Path baseDir

    Path homeDir

    Path currentDir

    boolean showAllProfiles

    String profile = DEFAULT_PROFILE

    boolean validateProfile

    List<Path> userConfigFiles = []

    List<Path> parsedConfigFiles = []

    List<String> parsedProfileNames

    boolean showClosures

    boolean showMissingVariables

    Map<ConfigObject, String> emptyVariables = new LinkedHashMap<>(10)

    Map<String,String> env = new HashMap<>(System.getenv())

    List<String> warnings = new ArrayList<>(10);

    {
        setHomeDir(Const.APP_HOME_DIR)
        setCurrentDir(Paths.get('.'))
    }

    ConfigBuilder setShowClosures(boolean value) {
        this.showClosures = value
        return this
    }

    ConfigBuilder showMissingVariables(boolean value) {
        this.showMissingVariables = value
        return this
    }

    ConfigBuilder setOptions( CliOptions options ) {
        this.options = options
        return this
    }

    ConfigBuilder setCmdRun( CmdRun cmdRun ) {
        this.cmdRun = cmdRun
        setProfile(cmdRun.profile)
        return this
    }

    ConfigBuilder setBaseDir( Path path ) {
        this.baseDir = path.complete()
        return this
    }

    ConfigBuilder setCurrentDir( Path path ) {
        this.currentDir = path.complete()
        return this
    }

    ConfigBuilder setHomeDir( Path path ) {
        this.homeDir = path.complete()
        return this
    }

    ConfigBuilder setCmdNode( CmdNode node ) {
        this.cmdNode = node
        return this
    }

    ConfigBuilder setCmdConfig( CmdConfig cmdConfig ) {
        showAllProfiles = cmdConfig.showAllProfiles
        setProfile(cmdConfig.profile)
        return this
    }

    ConfigBuilder setProfile( String value ) {
        profile = value ?: DEFAULT_PROFILE
        validateProfile = value as boolean
        return this
    }

    ConfigBuilder setShowAllProfiles(boolean value) {
        this.showAllProfiles = value
        return this
    }

    ConfigBuilder setUserConfigFiles( Path... files )  {
        setUserConfigFiles(files as List)
        return this
    }

    ConfigBuilder setUserConfigFiles( List<Path> files ) {
        if( files )
            userConfigFiles.addAll(files)
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
    @PackageScope
    List<Path> validateConfigFiles( List<String> files ) {

        def result = []
        if ( files ) {
            for( String fileName : files ) { 
                def thisFile = currentDir.resolve(fileName)
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
        def home = homeDir.resolve('config')
        if( home.exists() ) {
            log.debug "Found config home: $home"
            result << home
        }

        /**
         * Config file in the pipeline base dir
         * This config file name should be predictable, therefore cannot be overridden
         */
        def base = null
        if( baseDir && baseDir != currentDir ) {
            base = baseDir.resolve('nextflow.config')
            if( base.exists() ) {
                log.debug "Found config base: $base"
                result << base
            }
        }

        /**
         * Local or user provided file
         * Default config file name can be overridden with `NXF_CONFIG_FILE` env variable
         */
        def configFileName = env.get('NXF_CONFIG_FILE') ?: 'nextflow.config'
        def local = currentDir.resolve(configFileName)
        if( local.exists() && local != base ) {
            log.debug "Found config local: $local"
            result << local
        }

        def customConfigs = []
        if( userConfigFiles ) customConfigs.addAll(userConfigFiles)
        if( options?.userConfig ) customConfigs.addAll(options.userConfig)
        if( cmdRun?.runConfig ) customConfigs.addAll(cmdRun.runConfig)
        if( customConfigs ) {
            for( def item : customConfigs ) {
                def configFile = item instanceof Path ? item : currentDir.resolve(item.toString())
                if(!configFile.exists()) {
                    throw new AbortOperationException("The specified configuration file does not exist: $configFile -- check the name or choose another file")
                }

                log.debug "User config file: $configFile"
                result << configFile
            }
        }

        return result
    }

    /**
     * Create the nextflow configuration {@link ConfigObject} given a one or more
     * config files
     *
     * @param files A list of config files {@link Path}
     * @return The resulting {@link ConfigObject} instance
     */
    @PackageScope
    ConfigObject buildGivenFiles(List<Path> files) {

        final Map<String,String> vars = cmdRun?.env
        final boolean exportSysEnv = cmdRun?.exportSysEnv

        def items = []
        if( files ) for( Path file : files ) {
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

    @PackageScope
    ConfigObject buildGivenFiles(Path... files) {
        buildGivenFiles(files as List<Path>)
    }

    protected ConfigObject buildConfig0( Map env, List configEntries )  {
        assert env != null

        final slurper = new ConfigParser().setRenderClosureAsString(showClosures)
        ConfigObject result = new ConfigObject()

        if( cmdRun && (cmdRun.params || cmdRun.paramsFile) )
            slurper.setParams(cmdRun.parsedParams)

        // add the user specified environment to the session env
        env.sort().each { name, value -> result.env.put(name,value) }

        if( configEntries ) {
            // the configuration object binds always the current environment
            // so that in the configuration file may be referenced any variable
            // in the current environment
            final binding = new HashMap(System.getenv())
            binding.putAll(env)
            binding.put('baseDir', baseDir)
            binding.put('projectDir', baseDir)
            binding.put('launchDir', Paths.get('.').toRealPath())

            slurper.setBinding(binding)

            // merge of the provided configuration files
            for( def entry : configEntries ) {

                try {
                    merge0(result, slurper, entry)
                }
                catch( ConfigParseException e ) {
                    throw e
                }
                catch( Exception e ) {
                    def message = (entry instanceof Path ? "Unable to parse config file: '$entry'" : "Unable to parse configuration ")
                    throw new ConfigParseException(message,e)
                }
            }

            this.parsedProfileNames = new ArrayList<>(slurper.getProfileNames())
            if( validateProfile ) {
                checkValidProfile(slurper.getConditionalBlockNames())
            }

        }

        // guarantee top scopes
        for( String name : ['env','session','params','process','executor']) {
            if( !result.isSet(name) ) result.put(name, new ConfigObject())
        }

        return result
    }

    /**
     * Merge the main config with a separate config file
     *
     * @param result The main {@link ConfigObject}
     * @param slurper The {@ComposedConfigSlurper} parsed instance
     * @param entry The next config snippet/file to be parsed
     * @return
     */
    protected void merge0(ConfigObject result, ConfigParser slurper, entry) {
        if( !entry )
            return

        // select the profile
        if( showAllProfiles ) {
            def config = parse0(slurper,entry)
            validate(config,entry)
            result.merge(config)
            return
        }

        log.debug "Applying config profile: `${profile}`"
        def allNames = profile.tokenize(',')
        slurper.registerConditionalBlock('profiles', allNames)

        def config = parse0(slurper,entry)
        validate(config,entry)
        result.merge(config)
    }

    protected ConfigObject parse0(ConfigParser slurper, entry) {
        if( entry instanceof File ) {
            final path = entry.toPath()
            parsedConfigFiles << path
            return slurper.parse(path)
        }

        if( entry instanceof Path ) {
            parsedConfigFiles << entry
            return slurper.parse(entry)
        }

        return slurper.parse(entry.toString())
    }

    /**
     * Validate a config object verifying is does not contains unresolved attributes
     *
     * @param config The {@link ConfigObject} to verify
     * @param file The source config file/snippet
     * @return
     */
    protected validate(ConfigObject config, file, String parent=null, List stack = new ArrayList()) {
        for( String key : new ArrayList<>(config.keySet()) ) {
            final value = config.get(key)
            if( value instanceof ConfigObject ) {
                final fqKey = parent ? "${parent}.${key}": key as String
                if( value.isEmpty() ) {
                    final msg = "Unknown config attribute `$fqKey` -- check config file: $file".toString()
                    if( showMissingVariables ) {
                        emptyVariables.put(value, key)
                        warnings.add(msg)
                    }
                    else {
                        log.debug("In the following config snippet the attribute `$fqKey` is empty:\n${->config.prettyPrint().indent('  ')}")
                        throw new ConfigParseException(msg)
                    }
                }
                else {
                    stack.push(config)
                    try {
                        if( !stack.contains(value)) {
                            validate(value, file, fqKey, stack)
                        }
                        else {
                            log.debug("Found a recursive config property: `$fqKey`")
                        }
                    }
                    finally {
                        stack.pop()
                    }
                }
            }
            else if( value instanceof GString && showMissingVariables ) {
                final str = (GString) value
                for( int i=0; i<str.values.length; i++ ) {
                    // try replace empty interpolated strings with variable handle
                    final arg = str.values[i]
                    final name = emptyVariables.get(arg)
                    if( name )
                        str.values[i] = '$' + name
                }
            }
        }
    }

    protected void checkValidProfile(Collection<String> validNames) {
        if( !profile || profile == DEFAULT_PROFILE ) {
            return
        }

        log.debug "Available config profiles: $validNames"
        for( String name : profile.tokenize(',') ) {
            if( name in validNames )
                continue

            def message = "Unknown configuration profile: '${name}'"
            def choices = validNames.closest(name)
            if( choices ) {
                message += "\n\nDid you mean one of these?\n"
                choices.each { message += "    ${it}\n" }
                message += '\n'
            }

            throw new AbortOperationException(message)
        }
    }

    private String normalizeResumeId( String uniqueId ) {
        if( !uniqueId )
            return null

        if( uniqueId == 'last' || uniqueId == 'true' ) {
            uniqueId = HistoryFile.DEFAULT.getLast()?.sessionId

            if( !uniqueId ) {
                log.warn "It appears you have never run this project before -- Option `-resume` is ignored"
            }
        }

        return uniqueId
    }

    @PackageScope
    void configRunOptions(ConfigObject config, Map env, CmdRun cmdRun) {

        // -- set config options
        if( cmdRun.cacheable != null )
            config.cacheable = cmdRun.cacheable

        // -- set the run name
        if( cmdRun.runName )
            config.runName = cmdRun.runName

        if( cmdRun.stubRun )
            config.stubRun = cmdRun.stubRun

        // -- sets the working directory
        if( cmdRun.workDir )
            config.workDir = cmdRun.workDir

        else if( !config.workDir )
            config.workDir = env.get('NXF_WORK') ?: 'work'

        if( cmdRun.bucketDir )
            config.bucketDir = cmdRun.bucketDir

        // -- sets the library path
        if( cmdRun.libPath )
            config.libDir = cmdRun.libPath

        else if ( !config.isSet('libDir') && env.get('NXF_LIB') )
            config.libDir = env.get('NXF_LIB')

        // -- override 'process' parameters defined on the cmd line
        cmdRun.process.each { name, value ->
            config.process[name] = parseValue(value)
        }

        // -- apply the conda environment
        if( cmdRun.withConda ) {
            config.process.conda = cmdRun.withConda
        }

        // -- sets the resume option
        if( cmdRun.resume )
            config.resume = cmdRun.resume

        if( config.isSet('resume') )
            config.resume = normalizeResumeId(config.resume as String)

        // -- sets `dumpKeys` option
        if( cmdRun.dumpHashes )
            config.dumpHashes = cmdRun.dumpHashes

        if( cmdRun.dumpChannels )
            config.dumpChannels = cmdRun.dumpChannels.tokenize(',')

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

        // -- sets trace file options
        if( cmdRun.withTrace ) {
            if( !(config.trace instanceof Map) )
                config.trace = [:]
            config.trace.enabled = true
            if( cmdRun.withTrace != '-' )
                config.trace.file = cmdRun.withTrace
            else if( !config.trace.file )
                config.trace.file = TraceFileObserver.DEF_FILE_NAME
        }

        // -- sets report report options
        if( cmdRun.withReport ) {
            if( !(config.report instanceof Map) )
                config.report = [:]
            config.report.enabled = true
            if( cmdRun.withReport != '-' )
                config.report.file = cmdRun.withReport
            else if( !config.report.file )
                config.report.file = ReportObserver.DEF_FILE_NAME
        }

        // -- sets timeline report options
        if( cmdRun.withTimeline ) {
            if( !(config.timeline instanceof Map) )
                config.timeline = [:]
            config.timeline.enabled = true
            if( cmdRun.withTimeline != '-' )
                config.timeline.file = cmdRun.withTimeline
            else if( !config.timeline.file )
                config.timeline.file = TimelineObserver.DEF_FILE_NAME
        }

        // -- sets DAG report options
        if( cmdRun.withDag ) {
            if( !(config.dag instanceof Map) )
                config.dag = [:]
            config.dag.enabled = true
            if( cmdRun.withDag != '-' )
                config.dag.file = cmdRun.withDag
            else if( !config.dag.file )
                config.dag.file = GraphObserver.DEF_FILE_NAME
        }

        if( cmdRun.withNotification ) {
            if( !(config.notification instanceof Map) )
                config.notification = [:]
            if( cmdRun.withNotification in ['true','false']) {
                config.notification.enabled = cmdRun.withNotification == 'true'
            }
            else {
                config.notification.enabled = true
                config.notification.to = cmdRun.withNotification
            }
        }

        // -- sets the messages options
        if( cmdRun.withWebLog ) {
            if( !(config.weblog instanceof Map) )
                config.weblog = [:]
            config.weblog.enabled = true
            if( cmdRun.withWebLog != '-' )
                config.weblog.url = cmdRun.withWebLog
            else if( !config.weblog.url )
                config.weblog.url = WebLogObserver.DEF_URL
        }

        // -- sets tower options
        if( cmdRun.withTower ) {
            if( !(config.tower instanceof Map) )
                config.tower = [:]
            config.tower.enabled = true
            config.tower.endpoint = cmdRun.withTower
        }

        // -- add the command line parameters to the 'taskConfig' object
        if( cmdRun.params || cmdRun.paramsFile )
            config.params.putAll( cmdRun.parsedParams )

        if( cmdRun.withoutDocker && config.docker instanceof Map ) {
            // disable docker execution
            log.debug "Disabling execution in Docker contained as requested by cli option `-without-docker`"
            config.docker.enabled = false
        }

        if( cmdRun.withDocker ) {
            configContainer(config, 'docker', cmdRun.withDocker)
        }

        if( cmdRun.withPodman ) {
            configContainer(config, 'podman', cmdRun.withPodman)
        }

        if( cmdRun.withSingularity ) {
            configContainer(config, 'singularity', cmdRun.withSingularity)
        }

        if( cmdRun.withCharliecloud ) {
            configContainer(config, 'charliecloud', cmdRun.withCharliecloud)
        }
    }

    private void configContainer(ConfigObject config, String engine, def cli) {
        log.debug "Enabling execution in ${engine.capitalize()} container as requested by cli option `-with-$engine ${cmdRun.withDocker}`"

        if( !config.containsKey(engine) )
            config.put(engine, [:])

        if( !(config.get(engine) instanceof Map) )
            throw new AbortOperationException("Invalid `$engine` definition in the config file")

        def containerConfig = (Map)config.get(engine)
        containerConfig.enabled = true
        if( cli != '-' ) {
            // this is supposed to be a docker image name
            config.process.container = cli
        }
        else if( containerConfig.image ) {
            config.process.container = containerConfig.image
        }

        if( !hasContainerDirective(config.process) )
            throw new AbortOperationException("You have requested to run with ${engine.capitalize()} but no image were specified")

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

    ConfigObject buildConfigObject() {
        // -- configuration file(s)
        def configFiles = validateConfigFiles(options?.config)
        def config = buildGivenFiles(configFiles)

        if( cmdRun )
            configRunOptions(config, System.getenv(), cmdRun)

        return config
    }


    /**
     * @return A the application options hold in a {@code ConfigObject} instance
     */
    Map build() {
        buildConfigObject().toMap()
    }


}
