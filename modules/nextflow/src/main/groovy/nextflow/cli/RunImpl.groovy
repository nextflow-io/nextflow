/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.cli

import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.util.regex.Pattern

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import groovyx.gpars.GParsConfig
import nextflow.Const
import nextflow.NF
import nextflow.NextflowMeta
import nextflow.SysEnv
import nextflow.config.ConfigBuilder
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.script.ScriptFile
import nextflow.script.ScriptRunner
import nextflow.secret.SecretsLoader
import nextflow.util.CustomPoolFactory
import nextflow.util.HistoryFile
import org.yaml.snakeyaml.Yaml

/**
 * CLI `run` sub-command
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class RunImpl {

    interface Options extends IHubOptions {
        String getPipeline()
        List<String> getArgs()
        Map<String,String> getParams()

        String getBucketDir()
        Boolean getCacheable()
        Map<String,String> getClusterOptions()
        Integer getDeep()
        Boolean getDisableJobsCancellation()
        boolean getDsl1()
        boolean getDsl2()
        String getDumpChannels()
        boolean getDumpHashes()
        String getEntryName()
        Map<String,String> getEnv()
        Map<String,String> getExecutorOptions()
        boolean getExportSysEnv()
        boolean getLatest()
        String getLibPath()
        String getMainScript()
        boolean getOffline()
        String getParamsFile()
        String getPlugins()
        long getPollInterval()
        Integer getPoolSize()
        boolean getPreview()
        Map<String,String> getProcessOptions()
        String getProfile()
        Integer getQueueSize()
        String getResume()
        String getRevision()
        List<String> getRunConfig()
        String getRunName()
        boolean getStubRun()
        String getTest()
        String getWithApptainer()
        String getWithConda()
        Boolean getWithoutConda()
        String getWithCharliecloud()
        String getWithDag()
        String getWithDocker()
        boolean getWithoutDocker()
        String getWithFusion()
        boolean getWithMpi()
        String getWithNotification()
        String getWithPodman()
        boolean getWithoutPodman()
        String getWithReport()
        String getWithSingularity()
        String getWithSpack()
        Boolean getWithoutSpack()
        String getWithTimeline()
        String getWithTower()
        String getWithTrace()
        String getWithWave()
        String getWithWebLog()
        String getWorkDir()

        String getLauncherCliString()
        ILauncherOptions getLauncherOptions()

        void setRunName(String runName)
    }

    private Map<String,String> sysEnv = System.getenv()

    static {
        // install the custom pool factory for GPars threads
        GParsConfig.poolFactory = new CustomPoolFactory()
    }

    @Delegate
    private Options options

    RunImpl(Options options) {
        this.options = options
    }

    /* For testing purposes only */
    RunImpl() {}

    Boolean getDisableJobsCancellation() {
        options.disableJobsCancellation != null
            ? options.disableJobsCancellation
            : sysEnv.get('NXF_DISABLE_JOBS_CANCELLATION') as boolean
    }

    boolean getOffline() {
        options.offline || System.getenv('NXF_OFFLINE') as boolean
    }

    String getParamsFile() {
        options.paramsFile ?: sysEnv.get('NXF_PARAMS_FILE')
    }

    void run() {
        if( !pipeline )
            throw new AbortOperationException("No project name was specified")

        if( withPodman && withoutPodman )
            throw new AbortOperationException("Command line options `with-podman` and `without-podman` cannot be specified at the same time")

        if( withDocker && withoutDocker )
            throw new AbortOperationException("Command line options `with-docker` and `without-docker` cannot be specified at the same time")

        if( withConda && withoutConda )
            throw new AbortOperationException("Command line options `with-conda` and `without-conda` cannot be specified at the same time")

        if( withSpack && withoutSpack )
            throw new AbortOperationException("Command line options `with-spack` and `without-spack` cannot be specified at the same time")

        if( offline && latest )
            throw new AbortOperationException("Command line options `latest` and `offline` cannot be specified at the same time")

        if( dsl1 && dsl2 )
            throw new AbortOperationException("Command line options `dsl1` and `dsl2` cannot be specified at the same time")

        checkRunName()

        log.info "N E X T F L O W  ~  version ${Const.APP_VER}"
        Plugins.init()

        // -- specify the arguments
        final scriptFile = getScriptFile(pipeline)

        // create the config object
        final builder = new ConfigBuilder()
                .setLauncherOptions(launcherOptions)
                .setRunOptions(this)
                .setBaseDir(scriptFile.parent)
        final config = builder .build()

        // check DSL syntax in the config
        launchInfo(config, scriptFile)

        // check if NXF_ variables are set in nextflow.config
        checkConfigEnv(config)

        // -- load plugins
        final cfg = plugins ? [plugins: plugins.tokenize(',')] : config
        Plugins.load(cfg)

        // -- load secret provider
        if( SecretsLoader.isEnabled() ) {
            final provider = SecretsLoader.instance.load()
            config.withSecretProvider(provider)
        }

        // -- create a new runner instance
        final runner = new ScriptRunner(config)
        runner.setScript(scriptFile)
        runner.setPreview(this.preview)
        runner.session.profile = profile
        runner.session.commandLine = launcherCliString
        runner.session.ansiLog = launcherOptions.ansiLog
        runner.session.disableJobsCancellation = getDisableJobsCancellation()

        final isTowerEnabled = config.navigate('tower.enabled') as Boolean
        if( isTowerEnabled || log.isTraceEnabled() )
            runner.session.resolvedConfig = ConfigBuilder.resolveConfig(scriptFile.parent, this)
        // note config files are collected during the build process
        // this line should be after `ConfigBuilder#build`
        runner.session.configFiles = builder.parsedConfigFiles
        // set the commit id (if any)
        runner.session.commitId = scriptFile.commitId
        if( this.test ) {
            runner.test(this.test, args)
            return
        }

        def info = InfoImpl.status( log.isTraceEnabled() )
        log.debug( '\n'+info )

        // -- add this run to the local history
        runner.verifyAndTrackHistory(launcherCliString, runName)

        // -- run it!
        runner.execute(args, this.entryName)
    }

    protected void checkRunName() {
        if( runName == 'last' )
            throw new AbortOperationException("Not a valid run name: `last`")
        if( runName && !matchRunName(runName) )
            throw new AbortOperationException("Not a valid run name: `$runName` -- It must match the pattern $RUN_NAME_PATTERN")

        if( !runName ) {
            if( HistoryFile.disabled() )
                throw new AbortOperationException("Missing workflow run name")
            // -- make sure the generated name does not exist already
            runName = HistoryFile.DEFAULT.generateNextName()
        }

        else if( !HistoryFile.disabled() && HistoryFile.DEFAULT.checkExistsByName(runName) )
            throw new AbortOperationException("Run name `$runName` has been already used -- Specify a different one")
    }

    static final public Pattern RUN_NAME_PATTERN = Pattern.compile(/^[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}$/, Pattern.CASE_INSENSITIVE)

    static protected boolean matchRunName(String name) {
        RUN_NAME_PATTERN.matcher(name).matches()
    }

    protected ScriptFile getScriptFile(String pipelineName) {
        try {
            getScriptFile0(pipelineName)
        }
        catch (IllegalArgumentException | AbortOperationException e) {
            if( e.message.startsWith("Not a valid project name:") && !guessIsRepo(pipelineName)) {
                throw new AbortOperationException("Cannot find script file: $pipelineName")
            }
            else
                throw e
        }
    }

    protected ScriptFile getScriptFile0(String pipelineName) {
        assert pipelineName

        /*
         * read from the stdin
         */
        if( pipelineName == '-' ) {
            def file = tryReadFromStdin()
            if( !file )
                throw new AbortOperationException("Cannot access `stdin` stream")

            if( revision )
                throw new AbortOperationException("Revision option cannot be used when running a script from stdin")

            return new ScriptFile(file)
        }

        /*
         * look for a file with the specified pipeline name
         */
        def script = new File(pipelineName)
        if( script.isDirectory()  ) {
            script = mainScript ? new File(mainScript) : new AssetManager().setLocalPath(script).getMainScriptFile()
        }

        if( script.exists() ) {
            if( revision )
                throw new AbortOperationException("Revision option cannot be used when running a local script")
            return new ScriptFile(script)
        }

        /*
         * try to look for a pipeline in the repository
         */
        def manager = new AssetManager(pipelineName, this)
        def repo = manager.getProject()

        boolean checkForUpdate = true
        if( !manager.isRunnable() || latest ) {
            if( offline )
                throw new AbortOperationException("Unknown project `$repo` -- NOTE: automatic download from remote repositories is disabled")
            log.info "Pulling $repo ..."
            def result = manager.download(revision,deep)
            if( result )
                log.info " $result"
            checkForUpdate = false
        }
        // checkout requested revision
        try {
            manager.checkout(revision)
            manager.updateModules()
            final scriptFile = manager.getScriptFile(mainScript)
            if( checkForUpdate && !offline )
                manager.checkRemoteStatus(scriptFile.revisionInfo)
            // return the script file
            return scriptFile
        }
        catch( AbortOperationException e ) {
            throw e
        }
        catch( Exception e ) {
            throw new AbortOperationException("Unknown error accessing project `$repo` -- Repository may be corrupted: ${manager.localPath}", e)
        }

    }

    static protected File tryReadFromStdin() {
        if( !System.in.available() )
            return null

        getScriptFromStream(System.in)
    }

    static protected File getScriptFromStream( InputStream input, String name = 'nextflow' ) {
        input != null
        File result = File.createTempFile(name, null)
        result.deleteOnExit()
        input.withReader { Reader reader -> result << reader }
        return result
    }

    static protected boolean guessIsRepo(String name) {
        if( FileHelper.getUrlProtocol(name) != null )
            return true
        if( name.startsWith('/') )
            return false
        if( name.startsWith('./') || name.startsWith('../') )
            return false
        if( name.endsWith('.nf') )
            return false
        if( name.count('/') != 1 )
            return false
        return true
    }

    static final public DSL1 = '1'
    static final public DSL2 = '2'

    protected void launchInfo(ConfigMap config, ScriptFile scriptFile) {
        // -- determine strict mode
        final defStrict = sysEnv.get('NXF_ENABLE_STRICT') ?: false
        final strictMode = config.navigate('nextflow.enable.strict', defStrict)
        if( strictMode ) {
            log.debug "Enabling nextflow strict mode"
            NextflowMeta.instance.strictMode(true)
        }
        // -- determine dsl mode
        final dsl = detectDslMode(config, scriptFile.main.text, sysEnv)
        NextflowMeta.instance.enableDsl(dsl)
        // -- show launch info
        final ver = NF.dsl2 ? DSL2 : DSL1
        final repo = scriptFile.repository ?: scriptFile.source
        final head = preview ? "* PREVIEW * $scriptFile.repository" : "Launching `$repo`"
        if( scriptFile.repository )
            log.info "${head} [$runName] DSL${ver} - revision: ${scriptFile.revisionInfo}"
        else
            log.info "${head} [$runName] DSL${ver} - revision: ${scriptFile.getScriptId()?.substring(0,10)}"
    }

    static String detectDslMode(ConfigMap config, String scriptText, Map sysEnv) {
        // -- try determine DSL version from config file

        final dsl = config.navigate('nextflow.enable.dsl') as String

        // -- script can still override the DSL version
        final scriptDsl = NextflowMeta.checkDslMode(scriptText)
        if( scriptDsl ) {
            log.debug("Applied DSL=$scriptDsl from script declararion")
            return scriptDsl
        }
        else if( dsl ) {
            log.debug("Applied DSL=$dsl from config declaration")
            return dsl
        }
        // -- if still unknown try probing for DSL1
        if( NextflowMeta.probeDsl1(scriptText) ) {
            log.debug "Applied DSL=1 by probing script field"
            return DSL1
        }

        final envDsl = sysEnv.get('NXF_DEFAULT_DSL')
        if( envDsl ) {
            log.debug "Applied DSL=$envDsl from NXF_DEFAULT_DSL variable"
            return envDsl
        }
        else {
            log.debug "Applied DSL=2 by global default"
            return DSL2
        }
    }

    protected checkConfigEnv(ConfigMap config) {
        // Warn about setting NXF_ environment variables within env config scope
        final env = config.env as Map<String, String>
        for( String name : env.keySet() ) {
            if( name.startsWith('NXF_') && name!='NXF_DEBUG' ) {
                final msg = "Nextflow variables must be defined in the launching environment - The following variable set in the config file is going to be ignored: '$name'"
                log.warn(msg)
            }
        }
    }

    /**
     * Check whether pipeline params were provided via CLI options
     * or params file.
     */
    boolean hasParams() {
        return getParams() || getParamsFile()
    }

    /**
     * Collect all pipeline params provided via CLI options and params file.
     *
     * Memoized to avoid parsing the same file multiple times.
     *
     * @param configVars
     */
    @Memoized
    Map parsedParams(Map configVars) {

        final result = [:]
        final file = getParamsFile()
        if( file ) {
            def path = validateParamsFile(file)
            def type = path.extension.toLowerCase() ?: null
            if( type == 'json' )
                readJsonFile(path, configVars, result)
            else if( type == 'yml' || type == 'yaml' )
                readYamlFile(path, configVars, result)
        }

        // set the CLI params
        if( !params )
            return result

        for( Map.Entry<String,String> entry : params ) {
            addParam( result, entry.key, entry.value )
        }
        return result
    }

    static final public List<String> VALID_PARAMS_FILE = ['json', 'yml', 'yaml']

    private Path validateParamsFile(String file) {

        def result = FileHelper.asPath(file)
        def ext = result.getExtension()
        if( !VALID_PARAMS_FILE.contains(ext) )
            throw new AbortOperationException("Not a valid params file extension: $file -- It must be one of the following: ${VALID_PARAMS_FILE.join(',')}")

        return result
    }

    private void readJsonFile(Path file, Map configVars, Map result) {
        try {
            def text = configVars ? replaceVars0(file.text, configVars) : file.text
            def json = (Map)new JsonSlurper().parseText(text)
            result.putAll(json)
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exists: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()} - Cause: ${e.message}", e)
        }
    }

    private void readYamlFile(Path file, Map configVars, Map result) {
        try {
            def text = configVars ? replaceVars0(file.text, configVars) : file.text
            def yaml = (Map)new Yaml().load(text)
            result.putAll(yaml)
        }
        catch (NoSuchFileException | FileNotFoundException e) {
            throw new AbortOperationException("Specified params file does not exists: ${file.toUriString()}")
        }
        catch( Exception e ) {
            throw new AbortOperationException("Cannot parse params file: ${file.toUriString()}", e)
        }
    }

    static private Pattern PARAMS_VAR = ~/(?m)\$\{(\p{javaJavaIdentifierStart}\p{javaJavaIdentifierPart}*)}/

    static protected String replaceVars0(String content, Map binding) {
        content.replaceAll(PARAMS_VAR) { List<String> matcher ->
            // - the regex matcher is represented as list
            // - the first element is the matching string ie. `${something}`
            // - the second element is the group content ie. `something`
            // - make sure the regex contains at least a group otherwise the closure
            // parameter is a string instead of a list of the call fail
            final placeholder = matcher.get(0)
            final key = matcher.get(1)

            if( !binding.containsKey(key) ) {
                final msg = "Missing params file variable: $placeholder"
                if(NF.strictMode)
                    throw new AbortOperationException(msg)
                log.warn msg
                return placeholder
            }

            return binding.get(key)
        }
    }

    static final private Pattern DOT_ESCAPED = ~/\\\./
    static final private Pattern DOT_NOT_ESCAPED = ~/(?<!\\)\./

    static protected void addParam(Map params, String key, String value, List path=[], String fullKey=null) {
        if( !fullKey )
            fullKey = key
        final m = DOT_NOT_ESCAPED.matcher(key)
        if( m.find() ) {
            final p = m.start()
            final root = key.substring(0, p)
            if( !root ) throw new AbortOperationException("Invalid parameter name: $fullKey")
            path.add(root)
            def nested = params.get(root)
            if( nested == null ) {
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            else if( nested !instanceof Map ) {
                log.warn "Command line parameter --${path.join('.')} is overwritten by --${fullKey}"
                nested = new LinkedHashMap<>()
                params.put(root, nested)
            }
            addParam((Map)nested, key.substring(p+1), value, path, fullKey)
        }
        else {
            params.put(key.replaceAll(DOT_ESCAPED,'.'), parseParamValue(value))
        }
    }

    static protected parseParamValue(String str) {
        if ( SysEnv.get('NXF_DISABLE_PARAMS_TYPE_DETECTION') )
            return str

        if ( str == null ) return null

        if ( str.toLowerCase() == 'true') return Boolean.TRUE
        if ( str.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str==~/\d+(\.\d+)?/ && str.isInteger() ) return str.toInteger()
        if ( str==~/\d+(\.\d+)?/ && str.isLong() ) return str.toLong()
        if ( str==~/\d+(\.\d+)?/ && str.isDouble() ) return str.toDouble()

        return str
    }

}
