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
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Const
import nextflow.Nextflow
import nextflow.Session
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

    private Map params = [:]

    /**
     * Unnamed arguments passed to script into the 'args' script variable
     */
    private String[] args

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


    CliRunner setParam( String name, String value ) {
        assert name
        this.params.put(name, parseValue(value))
        return this
    }

    CliRunner setParam( Map<String,String> params ) {
        assert params != null
        params.each { name, value -> setParam(name, value) }
        return this
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

    def execute( File scriptFile, String... args ) {
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

    def execute( String scriptText, String... args ) {
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

    protected AbstractScript parseScript( File file, String... args) {
        assert file

        this.scriptFile = file
        parseScript( file.text, args )
    }

    protected AbstractScript parseScript( String scriptText, String... args = null) {

        bindings.setArgs( args )
        bindings.setParams( params )

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addStaticStars(Nextflow.getName())

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
        jcommander = new JCommander(result,args)

        return result
    }

    /**
     * Program entry method
     *
     * @param args The program options as specified by the user on the CLI
     */
    public static void main(String[] args)  {

        // parse the program arguments - and - configure the logger
        def options = parseMainArgs(args)
        LoggerHelper.configureLogger( options.quiet, options.debug, options.trace )

        // -- print out the version number, then exit
        if ( options.version ) {
            println getVersion(true)
            System.exit(0)
        }

        if ( !options.quiet ) {
            println Const.LOGO
        }

        // -- print out the program help, then exit
        if( options.help || !options.arguments ) {
            jcommander.usage()
            System.exit(0)
        }

        // -- check script name
        File scriptFile = new File(options.arguments[0])
        if ( !scriptFile.exists() ) {
            System.err.println "The specified script file does not exist: '$scriptFile'"
            System.exit(1)
        }


        try {
            // -- configuration file(s)
            def configFiles = validateConfigFiles(options.config)
            def config = buildConfig(configFiles)

            // -- create a new runner instance
            def runner = new CliRunner(config)
            runner.cacheable = options.cacheable

            // -- define the variable bindings
            if ( options.params ) {
                runner.setParam( options.params )
            }

            // -- specify the arguments
            def scriptArgs = options.arguments.size()>1 ? options.arguments[1..-1] : null

            // -- run it!
            runner.execute(scriptFile, scriptArgs as String[])
        }

        catch( Throwable fail ) {
            log.error fail.message, fail
            System.exit(1)
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

        ConfigObject result = new ConfigSlurper().parse('env{}')

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
                "${APP_NAME.capitalize()} ${APP_VER}_${APP_BUILDNUM} - Build at ${new Date(APP_TIMESTAMP).format(DATETIME_FORMAT)}"
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

}
