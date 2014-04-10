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

import java.lang.reflect.Field
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.spi.FileSystemProvider

import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import groovy.transform.InheritConstructors
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Const
import nextflow.ExitCode
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.daemon.DaemonLauncher
import nextflow.exception.ConfigParseException
import nextflow.exception.InvalidArgumentException
import nextflow.exception.MissingLibraryException
import nextflow.executor.ServiceName
import nextflow.util.FileHelper
import nextflow.util.HistoryFile
import nextflow.util.LoggerHelper
import nextflow.util.ServiceDiscover
import org.apache.commons.io.FilenameUtils
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer

/**
 * Application main class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CliRunner {

    @InheritConstructors
    static class CliArgumentException extends RuntimeException { }

    /**
     * The underlying execution session
     */
    private Session session

    /**
     * The interpreted script object
     */
    private BaseScript script

    /**
     * The extra binding variables specified by the user
     */
    private CliBinding bindings

    /**
     * The script file
     */
    private File scriptFile

    /*
     * the script raw output
     */
    private def output

    /**
     * The script result
     */
    private def result

    /**
     * Class path extension, they may be JAR files or directories containing java/groovy classes
     */
    private List<File> libraries

    /**
     * Instantiate the runner object creating a new session
     */
    def CliRunner( ) {
        this( [:] )
    }

    def CliRunner( Map config ) {
        session = new Session(config)
        bindings = new CliBinding(session)
    }

    /**
     * Only for testing purposes
     *
     * @param config
     */
    def CliRunner( ConfigObject config ) {
        this(configToMap(config))
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

    def void setLibPath( String str ) {
        if( !str ) { return }

        def files = str.split( File.pathSeparator ).collect { new File(it) }
        files?.each { File file ->
            if( !file.exists() ) { throw new MissingLibraryException("Cannot find specified library: ${file.absolutePath}")  }
            addLibPaths(file)
        }

    }

    /**
     * Add a library to the list of classpath extension. When {@code path} refers to
     * a 'directory', other than the directory, try to add all '.jar' files
     * contained in the directory itself.
     *
     * @param path
     */
    protected void addLibPaths( File path ) {
        assert path

        if( libraries == null ) { libraries = [] }

        if( path.isFile() && path.name.endsWith('.jar') ) {
            libraries << path
        }
        else if( path.isDirectory() ) {
            libraries << path
            path.eachFileMatch( ~/.+\.jar$/ ) { if(it.isFile()) this.libraries << it }
        }
    }


    /**
     * @return The interpreted script object
     */
    BaseScript getScript() { script }

    /**
     * @return The result produced by the script execution
     */
    def getResult() { result }

    /**
     * Execute a Nextflow script, it does the following:
     * <li>parse the script
     * <li>launch script execution
     * <li>await for all tasks completion
     *
     * @param scriptFile The file containing the script to be executed
     * @param args The arguments to be passed to the script
     * @return The result as returned by the {@code #run} method
     */

    def execute( File scriptFile, List<String> args = null ) {
        assert scriptFile

        // set the script name attribute
        session.scriptName = FilenameUtils.getBaseName(scriptFile.toString())

        // set the file name attribute
        this.scriptFile = scriptFile

        // get the script text and execute it
        execute( scriptFile.text, args )
    }


    /**
     * Execute a Nextflow script, it does the following:
     * <li>parse the script
     * <li>launch script execution
     * <li>await for all tasks completion
     *
     * @param scriptFile The file containing the script to be executed
     * @param args The arguments to be passed to the script
     * @return The result as returned by the {@code #run} method
     */

    def execute( String scriptText, List<String> args = null ) {
        assert scriptText

        try {
            // start session
            session.start()
            // parse the script
            script = parseScript(scriptText, args)
            // run the code
            run()
        }
        finally {
            terminate()
        }

        return result
    }

    /**
     * Test the name specified by the {@code methodName}
     *
     * @param scriptText
     * @param methodName
     * @param args
     */
    def test ( String scriptText, String methodName, List<String> args = null ) {
        assert scriptText
        assert methodName

        script = parseScript(scriptText, args)
        def values = args ? args.collect { parseValue(it) } : null

        def methodsToTest
        if ( methodName == '%all' ) {
            methodsToTest = script.metaClass.getMethods().findAll { it.name.startsWith('test') }.collect { it.name }.unique()
        }
        else {
            methodsToTest = methodName.split(',') *.trim()
        }

        if( !methodsToTest ) {
            log.info ('No test defined')
            return
        }

        for( String name: methodsToTest ) {
            try {
                def metaMethod = script.metaClass.getMethods().find { it.name == name }

                if( !metaMethod ) {
                    log.error "Unknown function '$name' -- Check the name spelling"
                }
                else {
                    def result = metaMethod.invoke(script, values as Object[] )
                    if( result != null ) {
                        log.info "SUCCESS: '$name' == $result"
                    }
                    else {
                        log.info "SUCCESS: '$name'"
                    }
                }
            }

            catch( AssertionError | Exception e ) {

                def lines = e.getMessage()?.readLines()?.collect { '  ' + it }
                def fmtMsg = lines ? "\n\n${lines.join('\n')}\n" : ''

                log.info "FAILED: function '$name'${fmtMsg}", e
            }
        }
    }

    def test( File scriptFile, String methodName, List<String> args = null ) {
        assert scriptFile
        assert methodName

        // set the script name attribute
        session.scriptName = FilenameUtils.getBaseName(scriptFile.toString())
        // set the file name attribute
        this.scriptFile = scriptFile

        test( scriptFile.text, methodName, args )

    }


    def normalizeOutput() {
        if( output instanceof Collection || output.getClass().isArray()) {
            result = (output as Collection)
        }
        else {
            result = output
        }
    }

    @Deprecated
    def normalizeOutput(def value) {
        if ( value instanceof DataflowReadChannel || value instanceof DataflowWriteChannel ) {
            return session.isAborted() ? null : Nextflow.read(value)
        }
        else {
            return value
        }
    }

    protected BaseScript parseScript( File file, List<String> args = null ) {
        assert file

        this.scriptFile = file
        parseScript( file.text, args )
    }

    protected BaseScript parseScript( String scriptText, List<String> args = null) {

        bindings.setArgs( new ArgsList(args) )
        bindings.setParams( session.config.params as Map )

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, groovy.transform.Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addStaticStars( Nextflow.name )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'splitter' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'operator' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'prioritySelector' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'select' )
        importCustomizer.addStaticImport( groovyx.gpars.dataflow.Dataflow.name, 'selector' )

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))

        // extend the class-loader if required
        def gcl = new GroovyClassLoader()
        if( libraries == null ) {
            // if no user defined libraries are provided
            // try to add the 'lib' folder in the local path if exist
            def localLib = session.baseDir ? new File(session.baseDir,'lib') : new File('lib')
            if( localLib.exists() ) {
                log.debug "Using default localLib path: $localLib"
                addLibPaths(localLib)
            }
        }

        libraries?.each { File lib -> def path = lib.absolutePath
            log.debug "Adding to the classpath library: ${path}"
            gcl.addClasspath(path)
            session.addClasspath(path)
        }

        // set the bytecode target directory
        def targetDir = FileHelper.createTempFolder(session.workDir).toFile()
        config.setTargetDirectory(targetDir)
        session.addClasspath(targetDir)

        // run and wait for termination
        BaseScript result
        def groovy = new GroovyShell(gcl, bindings, config)
        if ( scriptFile ) {
            result = groovy.parse( scriptText, scriptFile?.toString() ) as BaseScript
        }
        else {
            result = groovy.parse( scriptText ) as BaseScript
        }

        session.scriptClassName = result.class.name
        return result
    }


    /**
     * Launch the Nextflow script execution
     *
     * @return The value as returned by the user provided script
     */
    protected run() {
        assert script, "Missing script instance to run"

        // -- launch the script execution
        output = script.run()
    }

    protected terminate() {
        session.await()
        normalizeOutput()
        session.destroy()
    }


    /**
     * Create the application command line parser
     *
     * @return An instance of {@code CliBuilder}
     */

    static JCommander jcommander

    static boolean fullVersion

    @PackageScope
    static CliOptions parseMainArgs(String... args) {

        def options = CliOptions.normalizeArgs(args)
        def result = new CliOptions()
        jcommander = new JCommander(result, options as String[] )
        jcommander.setProgramName( Const.APP_NAME )
        fullVersion = '-version' in options
        return result
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
            def options = parseMainArgs(args)
            LoggerHelper.configureLogger(options)
            log.debug "${Const.APP_NAME} ${args?.join(' ')}"

            // -- print out the version number, then exit
            if ( options.version ) {
                println getVersion(fullVersion)
                System.exit(ExitCode.OK)
            }

            // -- print the history of executed commands
            if( options.history ) {
                HistoryFile.history.print()
                System.exit(ExitCode.OK)
            }

            // -- print out the program help, then exit
            if( options.help ) {
                println Const.LOGO
                jcommander.usage()
                System.exit(ExitCode.OK)
            }

            log.debug "${getVersion()} ~ host: ${localNameAndAddress}"
            File scriptFile = null
            if( !options.isDaemon() ) {
                // -- the script is the first item in the args
                if( options.arguments ) {
                    try {
                        scriptFile = getScriptFile(options.arguments[0])
                    }
                    catch( IOException e ) {
                        log.error "The specified script does not exist: '${options.arguments[0]}'\n", e
                        System.exit( ExitCode.MISSING_SCRIPT_FILE )
                    }

                    scriptBaseName = scriptFile.getBaseName()
                }
                // -- try to get the script from the stdin
                else {
                    scriptFile = tryReadFromStdin()
                    if( !scriptFile || scriptFile.empty() ) {
                        println Const.LOGO
                        jcommander.usage()
                        System.exit(ExitCode.OK)
                    }
                }

            }

            if( !options.quiet ) {
                println "N E X T F L O W  ~  version ${Const.APP_VER}"
            }

            // -- check file system providers
            checkFileSystemProviders()

            // create the config object
            def config = makeConfig(options)

            // -- launch daemon
            if( options.isDaemon() ) {
                log.debug "Launching cluster daemon"
                launchDaemon(config)
                return
            }

            // -- create a new runner instance
            def runner = new CliRunner(config)
            runner.session.cacheable = options.cacheable
            runner.session.resumeMode = options.resume != null
            runner.session.workDir = Paths.get(options.workDir).toAbsolutePath()
            runner.session.baseDir = scriptFile?.canonicalFile?.parentFile
            runner.libPath = options.libPath

            log.debug "Script bin dir: ${runner.session.binDir}"

            // -- specify the arguments
            def scriptArgs = options.arguments?.size()>1 ? options.arguments[1..-1] : []

            if( options.test ) {
                runner.test(scriptFile, options.test, scriptArgs )
            }
            else {
                // -- add this run to the local history
                HistoryFile.history.append( runner.session.uniqueId, args )
                // -- run it!
                runner.execute(scriptFile,scriptArgs)
            }
        }

        catch( ParameterException e ) {
            // print command line parsing errors
            // note: use system.err.println since if an exception is raised
            //       parsing the cli params the logging is not configured
            System.err.println "${e.getMessage()} -- Check the available command line parameters and syntax using '-h'"
            System.exit( ExitCode.INVALID_COMMAND_LINE_PARAMETER )
        }

        catch( ConfigParseException e )  {
            log.error "${e.message}\n${e.cause}\n\n"
            System.exit( ExitCode.INVALID_CONFIG)
        }

        catch ( MissingPropertyException e ) {
            log.error LoggerHelper.getErrorMessage(e, scriptBaseName), e
            System.exit( ExitCode.MISSING_PROPERTY )
        }

        catch( Throwable fail ) {
            log.error fail.message, fail
            System.exit( ExitCode.UNKNOWN_ERROR )
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
        input.withReader { reader -> result << reader }
        return result
    }

    static protected File getScriptFile( String urlOrPath ) {
        def lower = urlOrPath.toLowerCase()
        def isUrl = ['http','https','ftp'].any { lower.startsWith(it+'://') }

        if( isUrl ) {
            def url = new URL(urlOrPath)
            def fileName = new File(url.getPath()).getBaseName()
            return getScriptFromStream( url.newInputStream(), fileName )
        }

        def file = new File(urlOrPath)
        if( !file.exists() ) {
            throw new FileNotFoundException("File do not exist: $file")
        }
        return file
    }

    static private void checkFileSystemProviders() {

        // check if this class has been loaded
        boolean isInstalled = false
        FileSystemProvider.installedProviders().each {
            log.debug "Installed File System: '${it.scheme}' [${it.class.simpleName}]"
            if( it.scheme == 'dxfs' ) {
                isInstalled = true
            }
        }

        if( !isInstalled ) {
            // try to load DnaNexus file system provider dynamically
            Class provider
            try {
                provider = Class.forName('nextflow.fs.dx.DxFileSystemProvider')
            }
            catch( ClassNotFoundException e ) {
                log.debug "DxFileSystemProvider NOT available"
                return
            }

            // add it manually
            Field field = FileSystemProvider.class.getDeclaredField('installedProviders')
            field.setAccessible(true)
            List installedProviders = new ArrayList((List)field.get(null))
            installedProviders.add( provider.newInstance() )
            field.set(this, Collections.unmodifiableList(installedProviders))
            log.debug "Added 'DxFileSystemProvider' to list of installed providers [dxfs]"
        }

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
    def static List<File> validateConfigFiles( List<String> files ) {

        def result = []
        if ( files ) {
            files.each { String fileName ->
                def thisFile = new File(fileName)
                if(!thisFile.exists()) {
                    throw new CliArgumentException("The specified configuration file does not exist: $thisFile -- check the name or choose another file")
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


    def static Map buildConfig( List<File> files, CliOptions options ) {

        final Map<String,String> vars = options.env
        final boolean exportSysEnv = options.exportSysEnv

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

        if( options.executorOptions )  {
            def str = new StringBuilder()
            options.executorOptions.each { k, v ->
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


    static String getVersion(boolean full = false) {

        if ( full ) {
            Const.LOGO
        }
        else {
            "${Const.getAPP_NAME()} version ${Const.APP_VER}.${Const.APP_BUILDNUM}"
        }

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

    /**
     * Extends an {@code ArrayList} class adding a nicer index-out-of-range error message
     */
    static class ArgsList extends ArrayList<String> {


        ArgsList(List<String> values) {
            super( values ?: [] )
        }

        def String get( int pos ) {
            if( pos < 0 ) {
                throw new InvalidArgumentException("Argument array index cannot be lower than zero")
            }

            if( pos >= size() ) {
                throw new InvalidArgumentException("Arguments index out of range: $pos -- You may have not entered all arguments required by the pipeline")
            }

            super.get(pos)
        }

        def String getAt(int pos) {
            get(pos)
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
    static Map makeConfig( CliOptions options ) {

        // -- configuration file(s)
        def configFiles = validateConfigFiles(options.config)
        def config = buildConfig(configFiles, options)

        // -- override 'process' parameters defined on the cmd line
        options.process.each { name, value ->
            config.process[name] = parseValue(value)
        }

        // -- check for the 'continue' flag
        if( options.resume ) {
            def uniqueId = options.resume
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
        if( options.poolSize ) {
            config.poolSize = options.poolSize
        }
        if( options.queueSize ) {
            config.executor.queueSize = options.queueSize
        }
        if( options.pollInterval ) {
            config.executor.pollInterval = options.pollInterval
        }

        // -- add the command line parameters to the 'taskConfig' object
        options.params?.each { name, value ->
            config.params.put(name, parseValue(value))
        }

        return config
    }

    /**
     * Launch the daemon service
     *
     * @param config The nextflow configuration map
     */
    static launchDaemon( Map config ) {

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

    static private getLocalNameAndAddress() {
        def host = InetAddress.getLocalHost()
        "${host.getHostName()} [${host.getHostAddress()}]"
    }


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
            throw new IllegalArgumentException("Unknown daemon name: $name")

        clazz.newInstance()
    }

    static DaemonLauncher loadDaemonByClass( String name ) {
        try {
            return (DaemonLauncher)Class.forName(name).newInstance()
        }
        catch( Exception e ) {
            throw new IllegalStateException("Cannot load daemon: ${name}")
        }
    }

    static DaemonLauncher loadDaemonFirst() {
        def loader = ServiceLoader.load(DaemonLauncher).iterator()
        if( !loader.hasNext() )
            throw new IllegalStateException("No daemon services are available -- Cannot launch Nextflow in damon mode")

        return loader.next()
    }
}
