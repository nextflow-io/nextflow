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

package nextflow.script
import com.beust.jcommander.JCommander
import com.beust.jcommander.ParameterException
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Const
import nextflow.ExitCode
import nextflow.Nextflow
import nextflow.Session
import nextflow.exception.InvalidArgumentException
import nextflow.util.HistoryFile
import nextflow.util.LoggerHelper
import org.apache.commons.io.FilenameUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CliRunner {

    @InheritConstructors
    static class CliArgumentException extends RuntimeException {

    }

    /**
     * The underlying execution session
     */
    private Session session

    /**
     * The interpreted script object
     */
    private AbstractScript script

    /**
     * The directory where the output folder is created, optional.
     * If it is not specified, the output is stored in folder created in the current directory
     */
    private File workDirectory

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


    static private parseValue( String str ) {

        if ( str == null ) return null

        if ( str?.toLowerCase() == 'true') return Boolean.TRUE
        if ( str?.toLowerCase() == 'false' ) return Boolean.FALSE

        if ( str.isInteger() ) return str.toInteger()
        if ( str.isLong() ) return str.toLong()
        if ( str.isDouble() ) return str.toDouble()

        return str

    }

    /**
     * Enable/Disable tasks results caching for the current session
     *
     * @param value
     * @return
     */
    CliRunner setCacheable( boolean value ) {
        session.cacheable = value
        return this
    }

    /**
     * @return The interpreted script object
     */
    AbstractScript getScript() { script }

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

        script = parseScript(scriptText, args)
        try {
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
        if ( output == null ) {
            result = null
        }
        else if( output instanceof Collection || output.getClass().isArray()) {
            result = (output as Collection).collect { normalizeOutput(it) }
        }
        else {
            result = normalizeOutput(output)
        }
    }

    def normalizeOutput(def value) {
        if ( value instanceof DataflowReadChannel || value instanceof DataflowWriteChannel ) {
            return session.isAborted() ? null : Nextflow.read(value)
        }
        else {
            return value
        }
    }

    protected AbstractScript parseScript( File file, List<String> args = null ) {
        assert file

        this.scriptFile = file
        parseScript( file.text, args )
    }

    protected AbstractScript parseScript( String scriptText, List<String> args = null) {

        bindings.setArgs( new ArgsList(args) )
        bindings.setParams( session.config.params as Map )

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addStaticStars( Nextflow.getName() )
        importCustomizer.addImports( Channel.getName() )

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = AbstractScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskScriptClosureTransform))

        // run and wait for termination
        def groovy = new GroovyShell(this.class.classLoader, bindings, config)
        if ( scriptFile ) {
            groovy.parse( scriptText, scriptFile?.toString() ) as AbstractScript
        }
        else {
            groovy.parse( scriptText ) as AbstractScript
        }
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
        session?.terminate()
    }


    /**
     * Create the application command line parser
     *
     * @return An instance of {@code CliBuilder}
     */

    static JCommander jcommander

    static private CliOptions parseMainArgs(String... args) {

        def result = new CliOptions()
        jcommander = new JCommander(result, normalizeArgs( args ) as String[] )

        return result
    }

    static protected List<String> normalizeArgs( String ... args ) {

        def normalized = []
        int i=0
        while( true ) {
            if( i==args.size() ) { break }

            def current = args[i++]
            normalized << current


            if( current == '-resume' ) {
                if( i<args.size() && !args[i].startsWith('-') && (args[i]=='last' || args[i] =~~ /[0-9a-f]{8}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{4}\-[0-9a-f]{8}/) ) {
                    normalized << args[i++]
                }
                else {
                    normalized << 'last'
                }
            }
            else if( current == '-test' && (i==args.size() || args[i].startsWith('-'))) {
                normalized << '%all'
            }
        }

        return normalized
    }

    /**
     * Program entry method
     *
     * @param args The program options as specified by the user on the CLI
     */
    public static void main(String... args)  {

        try {
            // -- parse the program arguments - and - configure the logger
            def options = parseMainArgs(args)
            LoggerHelper.configureLogger( options.quiet, options.debug, options.trace )

            // -- print out the version number, then exit
            if ( options.version ) {
                println getVersion(true)
                System.exit(ExitCode.OK)
            }

            // -- print the history of executed commands
            if( options.history ) {
                HistoryFile.history.print()
                System.exit(ExitCode.OK)
            }

            if ( !options.quiet ) {
                println Const.LOGO
            }

            // -- print out the program help, then exit
            if( options.help || !options.arguments ) {
                jcommander.usage()
                System.exit(ExitCode.OK)
            }

            // -- check script name
            File scriptFile = new File(options.arguments[0])
            if ( !scriptFile.exists() ) {
                log.error "The specified script file does not exist: '$scriptFile'"
                System.exit( ExitCode.MISSING_SCRIPT_FILE )
            }

            // -- configuration file(s)
            def configFiles = validateConfigFiles(options.config)
            def config = buildConfig(configFiles)

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

            // -- add the command line parameters to the 'config' object
            options.params?.each { name, value ->
                config.params.put(name, parseValue(value))
            }

            // -- create a new runner instance
            def runner = new CliRunner(config)
            runner.cacheable = options.cacheable

            // -- specify the arguments
            def scriptArgs = options.arguments.size()>1 ? options.arguments[1..-1] : null

            if( options.test ) {
                runner.test(scriptFile, options.test, scriptArgs )
            }
            else {
                // -- set a shutdown hook to save the current session ID and command lines
                addShutdownHook { HistoryFile.history.append( runner.session.uniqueId, args ) }
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

        catch( Throwable fail ) {
            log.error fail.message, fail
            System.exit( ExitCode.UNKNOWN_ERROR )
        }

    }


    /**
     * Transform the specified list of string to a list of files, verifying their existence.
     * <p>
     *     If a file in the list does not exist an exception of type {@code CliArgumentException} is thrown.
     * <p>
     *     If the specified list is empty it tries to return of default configuration files located at:
     *     <li>$HOME/.nextflow/config
     *     <li>$PWD/nextflow.config
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


    def static Map buildConfig( List<File> files ) {

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

        buildConfig0( System.getenv(), texts )
    }


    def static Map buildConfig0( Map env, List<String> confText )  {
        assert env

        ConfigObject result = new ConfigSlurper().parse('env{}; session{}; params{} ')

        env.sort().each { name, value -> result.env.put(name,value) }

        confText?.each { String text ->

            if ( text ) {
                def cfg = new ConfigSlurper()
                cfg.setBinding(env)
                result.merge( cfg.parse(text) )
            }

        }

        // convert the ConfigObject to plain map
        // this because when accessing a non-existing entry in a ConfigObject it return and empty map as value
        return configToMap( result )
    }


    static String getVersion(boolean full = false) {

        Const.with {
            if ( full ) {
                "${getAPP_NAME()} version ${APP_VER}_${APP_BUILDNUM} - build timestamp ${new Date(APP_TIMESTAMP).format(DATETIME_FORMAT)}"
            }
            else {
                APP_VER
            }
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
}
