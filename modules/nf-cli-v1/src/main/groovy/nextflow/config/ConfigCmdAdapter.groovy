/*
 * Copyright 2013-2024, Seqera Labs
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

import java.nio.file.Path

import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.NF
import nextflow.cli.CliOptions
import nextflow.cli.CmdConfig
import nextflow.cli.CmdNode
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.trace.GraphObserver
import nextflow.trace.ReportObserver
import nextflow.trace.TimelineObserver
import nextflow.trace.TraceFileObserver
import nextflow.util.HistoryFile
import nextflow.util.SecretHelper

import static nextflow.util.ConfigHelper.toCanonicalString
/**
 * Adapter for building a Nextflow config with CLI overrides
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ConfigCmdAdapter {

    ConfigBuilder builder

    CliOptions options

    CmdRun cmdRun

    CmdNode cmdNode

    Path baseDir

    Path currentDir

    Path homeDir

    List<Path> userConfigFiles = []

    Map<String,String> sysEnv = new HashMap<>(System.getenv())

    ConfigCmdAdapter(ConfigBuilder builder) {
        this.builder = builder
        setHomeDir(Const.APP_HOME_DIR)
        setCurrentDir(Path.of('.'))
    }

    ConfigCmdAdapter() {
        this(new ConfigBuilder())
    }

    ConfigCmdAdapter setOptions(CliOptions options) {
        this.options = options
        builder.setIgnoreIncludes(options?.ignoreConfigIncludes)
        return this
    }

    ConfigCmdAdapter setCmdConfig(CmdConfig cmdConfig) {
        builder.setShowAllProfiles(cmdConfig.showAllProfiles)
        builder.setProfile(cmdConfig.profile)
        return this
    }

    ConfigCmdAdapter setCmdNode(CmdNode cmdNode) {
        this.cmdNode = cmdNode
        return this
    }

    ConfigCmdAdapter setCmdRun(CmdRun cmdRun) {
        this.cmdRun = cmdRun
        builder.setProfile(cmdRun.profile)
        return this
    }

    ConfigCmdAdapter setBaseDir(Path path) {
        this.baseDir = path.complete()
        builder.setBaseDir(baseDir)
        return this
    }

    ConfigCmdAdapter setCurrentDir(Path path) {
        this.currentDir = path.complete()
        builder.setCurrentDir(currentDir)
        return this
    }

    ConfigCmdAdapter setHomeDir(Path path) {
        this.homeDir = path.complete()
        return this
    }

    ConfigCmdAdapter setUserConfigFiles(Path... files)  {
        setUserConfigFiles(files as List<Path>)
        return this
    }

    ConfigCmdAdapter setUserConfigFiles(List<Path> files) {
        userConfigFiles.addAll(files)
        return this
    }

    /**
     * Build the configuration as a ConfigMap.
     */
    ConfigMap build() {
        toConfigMap(buildConfigObject())
    }

    /**
     * Build the configuration as a ConfigObject.
     */
    ConfigObject buildConfigObject() {
        // -- configuration file(s)
        final configFiles = resolveConfigFiles(options?.config)
        validateConfigFiles(configFiles)
        final config = buildGivenFiles(configFiles)

        if( cmdRun )
            configRunOptions(config, System.getenv(), cmdRun)

        return config
    }

    /**
     * Transform the specified list of string to a list of files.
     *
     * The following default configuration files are used if no files are specified:
     * - $HOME/.nextflow/config
     * - $PWD/nextflow.config
     *
     * @param files
     * @return
     */
    protected List<Path> resolveConfigFiles(List<String> files) {

        if( files ) {
            return files.stream()
                .map(filename -> currentDir.resolve(filename))
                .toList()
        }

        final List<Path> result = []

        /*
         * config file in the nextflow home
         */
        final home = homeDir.resolve('config')
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
        final configFileName = sysEnv.get('NXF_CONFIG_FILE') ?: 'nextflow.config'
        final local = currentDir.resolve(configFileName)
        if( local.exists() && local != base ) {
            log.debug "Found config local: $local"
            result << local
        }

        final customConfigs = []
        if( userConfigFiles )
            customConfigs.addAll(userConfigFiles)
        if( options?.userConfig )
            customConfigs.addAll(options.userConfig)
        if( cmdRun?.runConfig )
            customConfigs.addAll(cmdRun.runConfig)
        for( final item : customConfigs ) {
            final configFile = item instanceof Path ? item : currentDir.resolve(item.toString())
            log.debug "User config file: $configFile"
            result << configFile
        }

        return result
    }

    /**
     * Validate the existence of a list of config files.
     *
     * @param files
     */
    protected void validateConfigFiles(List<Path> files) {
        for( final file : files ) { 
            if( !file.exists() )
                throw new AbortOperationException("The specified configuration file does not exist: $file -- check the name or choose another file")
            if( !file.isFile() )
                throw new AbortOperationException("The specified configuration file is not a file: $file -- check the name or choose another file")
        }
    }

    /**
     * Create the nextflow configuration {@link ConfigObject} given a one or more
     * config files
     *
     * @param files A list of config files {@link Path}
     * @return The resulting {@link ConfigObject} instance
     */
    protected ConfigObject buildGivenFiles(List<Path> files) {

        if( cmdRun && (cmdRun.hasParams()) )
            builder.setParams(cmdRun.parsedParams(builder.configVars()))

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

        final Map<String,String> env = [:]
        if( cmdRun?.exportSysEnv ) {
            log.debug "Adding current system environment to session environment"
            env.putAll(System.getenv())
        }

        final envVars = cmdRun?.env
        if( envVars ) {
            log.debug "Adding the following variables to session environment: $envVars"
            env.putAll(envVars)
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

        return builder.build(env, items)
    }

    protected ConfigObject buildGivenFiles(Path... files) {
        return buildGivenFiles(files as List<Path>)
    }

    private static String wrapValue( value ) {
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
     * Apply command-line arguments to the config.
     *
     * @param config
     * @param env
     * @param cmdRun
     */
    @CompileDynamic
    protected void configRunOptions(ConfigObject config, Map env, CmdRun cmdRun) {

        // -- set config options
        if( cmdRun.cacheable != null )
            config.cacheable = cmdRun.cacheable

        // -- set the run name
        if( cmdRun.runName )
            config.runName = cmdRun.runName

        if( cmdRun.stubRun )
            config.stubRun = cmdRun.stubRun

        // -- set the output directory
        if( cmdRun.outputDir )
            config.outputDir = cmdRun.outputDir

        if( cmdRun.preview )
            config.preview = cmdRun.preview

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

        if( cmdRun.withoutConda && config.conda instanceof Map ) {
            // disable conda execution
            log.debug "Disabling execution with Conda as requested by command-line option `-without-conda`"
            config.conda.enabled = false
        }

        // -- apply the conda environment
        if( cmdRun.withConda ) {
            if( cmdRun.withConda != '-' )
                config.process.conda = cmdRun.withConda
            config.conda.enabled = true
        }

        if( cmdRun.withoutSpack && config.spack instanceof Map ) {
            // disable spack execution
            log.debug "Disabling execution with Spack as requested by command-line option `-without-spack`"
            config.spack.enabled = false
        }

        // -- apply the spack environment
        if( cmdRun.withSpack ) {
            if( cmdRun.withSpack != '-' )
                config.process.spack = cmdRun.withSpack
            config.spack.enabled = true
        }

        // -- sets the resume option
        if( cmdRun.resume )
            config.resume = cmdRun.resume

        if( config.isSet('resume') )
            config.resume = normalizeResumeId(config.resume as String)

        // -- sets `dumpHashes` option
        if( cmdRun.dumpHashes ) {
            config.dumpHashes = cmdRun.dumpHashes != '-' ? cmdRun.dumpHashes : 'default'
        }

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
            log.warn "The command line option '-with-weblog' is deprecated - consider enabling this feature by setting 'weblog.enabled=true' in your configuration file"
            if( !(config.weblog instanceof Map) )
                config.weblog = [:]
            config.weblog.enabled = true
            if( cmdRun.withWebLog != '-' )
                config.weblog.url = cmdRun.withWebLog
            else if( !config.weblog.url )
                config.weblog.url = 'http://localhost'
        }

        // -- sets tower options
        if( cmdRun.withTower ) {
            if( !(config.tower instanceof Map) )
                config.tower = [:]
            config.tower.enabled = true
            if( cmdRun.withTower != '-' )
                config.tower.endpoint = cmdRun.withTower
            else if( !config.tower.endpoint )
                config.tower.endpoint = 'https://api.cloud.seqera.io'
        }

        // -- set wave options
        if( cmdRun.withWave ) {
            if( !(config.wave instanceof Map) )
                config.wave = [:]
            config.wave.enabled = true
            if( cmdRun.withWave != '-' )
                config.wave.endpoint = cmdRun.withWave
            else if( !config.wave.endpoint )
                config.wave.endpoint = 'https://wave.seqera.io'
        }

        // -- set fusion options
        if( cmdRun.withFusion ) {
            if( !(config.fusion instanceof Map) )
                config.fusion = [:]
            config.fusion.enabled = cmdRun.withFusion == 'true'
        }

        // -- set cloudcache options
        final envCloudPath = env.get('NXF_CLOUDCACHE_PATH')
        if( cmdRun.cloudCachePath || envCloudPath ) {
            if( !(config.cloudcache instanceof Map) )
                config.cloudcache = [:]
            if( !config.cloudcache.isSet('enabled') )
                config.cloudcache.enabled = true
            if( cmdRun.cloudCachePath && cmdRun.cloudCachePath != '-' )
                config.cloudcache.path = cmdRun.cloudCachePath
            else if( !config.cloudcache.isSet('path') && envCloudPath )
                config.cloudcache.path = envCloudPath
        }

        // -- add the command line parameters to the config
        if( cmdRun.hasParams() )
            config.params = mergeMaps( (Map)config.params, cmdRun.parsedParams(builder.configVars()), NF.strictMode )

        if( cmdRun.withoutDocker && config.docker instanceof Map ) {
            // disable docker execution
            log.debug "Disabling execution in Docker container as requested by command-line option `-without-docker`"
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

        if( cmdRun.withApptainer ) {
            configContainer(config, 'apptainer', cmdRun.withApptainer)
        }

        if( cmdRun.withCharliecloud ) {
            configContainer(config, 'charliecloud', cmdRun.withCharliecloud)
        }
    }

    private String normalizeResumeId( String uniqueId ) {
        if( !uniqueId )
            return null
        if( uniqueId == 'last' || uniqueId == 'true' ) {
            if( HistoryFile.disabled() )
                throw new AbortOperationException("The resume session id should be specified via `-resume` option when history file tracking is disabled")
            uniqueId = HistoryFile.DEFAULT.getLast()?.sessionId

            if( !uniqueId ) {
                log.warn "It appears you have never run this project before -- Option `-resume` is ignored"
            }
        }

        return uniqueId
    }

    @CompileDynamic
    private void configContainer(ConfigObject config, String engine, def cli) {
        log.debug "Enabling execution in ${engine.capitalize()} container as requested by command-line option `-with-$engine ${cmdRun.withDocker}`"

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
            throw new AbortOperationException("You have requested to run with ${engine.capitalize()} but no image was specified")

    }

    /**
     * Verify that configuration for process contains at last one `container` directive
     *
     * @param process
     * @return {@code true} when a `container` is defined or {@code false} otherwise
     */
    @CompileDynamic
    protected boolean hasContainerDirective(process)  {

        if( process instanceof Map ) {
            if( process.container )
                return true

            def result = process
                    .findAll { String name, value -> (name.startsWith('withName:') || name.startsWith('$')) && value instanceof Map }
                    .find { String name, Map value -> value.container as boolean }  // the first non-empty `container` string

            return result as boolean
        }

        return false
    }

    /**
     * Merge two maps recursively avoiding keys to be overwritten
     *
     * @param config
     * @param params
     * @return a map resulting of merging result and right maps
     */
    protected Map mergeMaps(Map config, Map params, boolean strict, List keys=[]) {
        if( config==null )
            config = new LinkedHashMap()

        for( Map.Entry entry : params ) {
            final key = entry.key.toString()
            final value = entry.value
            final previous = getConfigVal0(config, key)
            keys << entry.key
            
            if( previous==null ) {
                config[key] = value
            }
            else if( previous instanceof Map && value instanceof Map ) {
                mergeMaps(previous, value, strict, keys)
            }
            else {
                if( previous instanceof Map || value instanceof Map ) {
                    final msg = "Configuration setting type with key '${keys.join('.')}' does not match the parameter with the same key - Config value=$previous; parameter value=$value"
                    if(strict)
                        throw new AbortOperationException(msg)
                    log.warn(msg)
                }
                config[key] = value
            }
        }

        return config
    }

    private Object getConfigVal0(Map config, String key) {
        if( config instanceof ConfigObject ) {
            return config.isSet(key) ? config.get(key) : null
        }
        else {
            return config.get(key)
        }
    }

    static String resolveConfig(Path baseDir, CmdRun cmdRun) {

        final builder = new ConfigBuilder()
                .setShowClosures(true)
                .setStripSecrets(true)

        final config = new ConfigCmdAdapter(builder)
                .setOptions(cmdRun.launcher.options)
                .setCmdRun(cmdRun)
                .setBaseDir(baseDir)
                .buildConfigObject()

        // strip secret
        SecretHelper.hideSecrets(config)
        // compute config
        final result = toCanonicalString(config, false)
        // dump config for debugging
        log.trace "Resolved config:\n${result.indent('\t')}"
        return result
    }

    protected static ConfigMap toConfigMap(ConfigObject config) {
        assert config != null
        (ConfigMap)normalize0((Map)config)
    }

    private static Object normalize0( config ) {

        if( config instanceof Map ) {
            ConfigMap result = new ConfigMap(config.size())
            for( String name : config.keySet() ) {
                def value = (config as Map).get(name)
                result.put(name, normalize0(value))
            }
            return result
        }
        else if( config instanceof Collection ) {
            List result = new ArrayList(config.size())
            for( entry in config ) {
                result << normalize0(entry)
            }
            return result
        }
        else {
            return config
        }
    }
}
