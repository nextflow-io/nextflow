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
import nextflow.Const
import nextflow.Nextflow
import nextflow.Session
import nextflow.util.FileHelper
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
    private CliBinding params = new CliBinding()

    /**
     * Unnamed arguments passed to script into the 'args' script variable
     */
    private String[] args

    /**
     * The script name (without the file extension)
     */
    private String name

    /**
     * The script file
     */
    private File scriptFile

    /**
     * Instantiate the runner object creating a new session
     */
    def CliRunner( ) {
        this( [:] )
    }

    def CliRunner( Map config ) {
        session = new Session(config)
    }

    /**
     * Only for testing purposes
     *
     * @param config
     */
    def CliRunner( ConfigObject config ) {
        this(configToMap(config))
    }


    /**
     * Instantiate the runner object using the specified session
     *
     * @param session
     */
    def CliRunner( Session session ) {
        this.session = session
    }


    CliRunner setParam( String name, String value ) {
        assert name
        this.params.setParam(name, value)
        return this
    }

    CliRunner setParam( Map<String,String> params ) {
        assert params != null
        params.each { name, value -> setParam(name, value) }
        return this
    }

    /**
     * @return The interpreted script object
     */
    AbstractScript getScript() { script }

    File getWorkDirectory() { workDirectory }

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
        if( !this.name ) {
            this.name = FilenameUtils.getBaseName(scriptFile.toString())
        }

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
            return run()
        }
        finally {
            terminate()
        }

    }


    protected AbstractScript parseScript( File file, String... args) {
        assert file

        this.scriptFile = file
        parseScript( file.text, args )
    }

    protected AbstractScript parseScript( String scriptText, String... args = null) {

        params['__$session'] = session
        params['args'] = args

//        // comment
//        def lines = scriptText.readLines()
//        if ( lines || lines[0].startsWith('#!') ) {
//            lines[0] = "//" + lines[0]
//        }

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addStaticStars(Nextflow.getName())

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = AbstractScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskScriptClosureTransform))

        // run and wait for termination
        def groovy = new GroovyShell(this.class.classLoader, params, config)
        if ( scriptFile ) {
            groovy.parse( scriptText, scriptFile?.toString() ) as AbstractScript
        }
        else {
            groovy.parse( scriptText ) as AbstractScript
        }
    }

    protected void defineWorkDirectory() {

        // -- the user may have specified teh work dir,
        //    if no create a new folder in the current directory
        if ( !workDirectory ) {
            workDirectory = FileHelper.tryCreateDir(new File('.', name ?: 'run'))
        }

        // -- make sure it is empty
        if ( FileHelper.isNotEmpty(workDirectory) ) {
            throw new IllegalStateException("Working directory must be empty: '${workDirectory.absolutePath}' -- Empty it or specify a different working directory to continue")
        }

        // -- se the 'run' name if missing
        if ( !name ) {
            name = FilenameUtils.getBaseName( workDirectory.name )
        }

        // -- set the working directory into the session
        session.workDirectory = workDirectory

    }


    /**
     * Launch the Nextflow script execution
     *
     * @return The value as returned by the user provided script
     */
    protected run() {
        assert script, "Missing script instance to run"

        // -- define work-directory
        defineWorkDirectory()

        // -- launch the script execution
        def result = script.run()

        // -- normalize the
        if ( result == null ) {
            return null
        }
        else if( result instanceof Collection || result.getClass().isArray()) {
            return result.collect { Nextflow.read(it) }
        }
        else {
            Nextflow.read(result)
        }


    }

    protected terminate() {
        session?.awaitTermination()
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
        LoggerHelper.configureLogger( options.debug, options.trace )

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
        String scriptName = FilenameUtils.getBaseName(options.arguments[0])
        if ( !scriptFile.exists() ) {
            System.err.println "The specified script file does not exist: '$scriptFile'"
            System.exit(1)
        }


        try {
            // -- configuration file(s)
            def configFiles = validateConfigFiles(options.config)
            def config = buildConfig(options.exportEnvironment, configFiles)

            // -- create a new runner instance
            def runner = new CliRunner(config)

            // -- define the working directory
            runner.workDirectory = validateWorkDirectory( options.workDirectory, scriptName )
            log.info "Work directory: ${runner.workDirectory}"

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
     * If the workDirectory string is specified, create that folder
     * otherwise tries to create a temporary working directory in current folder
     * using the 'scriptName' argument
     *
     * @param workDirectory The directory path that need to be created or {@code null} if a temporary working path is required
     * @param folderName The name of a temporary folder, used only when @{code workDirectory} is {@code null}
     * @return The directory created
     */
    static File validateWorkDirectory(String workDirectory, String folderName) {

        File workPath = null

        /*
         * verify user specified 'work directory'
         */
        if ( workDirectory ) {
            workPath = new File(workDirectory).canonicalFile

            if ( workPath.exists() ) {
                if ( !workPath.isDirectory() ) {
                    throw new CliArgumentException("The specified work directory is not valid: '${workDirectory}' -- Specify a valid path directory path")
                }

                if ( !FileHelper.isEmpty(workPath) ) {
                    throw new CliArgumentException("The work directory must be empty: '${workPath}' -- Remove the content or choose a different path")
                }
            }
            else {
                if ( !workPath.mkdirs() ) {
                    throw new CliArgumentException("Cannot create work directory: ${workPath} -- Check access permission or choose a different path")
                }
            }

            return workPath
        }

        /*
         * create a new working directory in the current folder using the script name as sub-folder name
         */
        FileHelper.tryCreateDir( new File('.',"run-${folderName}") )
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


    def static Map buildConfig( boolean export, List<File> files ) {

        def texts = []
        files?.each { File file ->
            if (!file.exists()) {
                log.warn "The specified configuration file cannot be found: $file"
            }
            else {
                texts << file.text
            }
        }

        buildConfig0( System.getenv(), export, texts )
    }


    def static Map buildConfig0( Map env, boolean export, List<String> confText )  {
        assert env

        ConfigObject result = new ConfigSlurper().parse('env{}')

        if ( export ) {
            env.sort().each { name, value -> result.env.put(name,value) }
        }

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
                "${APP_NAME.capitalize()} ${APP_VER}_${APP_BUILDNUM} - Build at ${new Date(Const.APP_TIMESTAMP).format(DATETIME_FORMAT)}"
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
