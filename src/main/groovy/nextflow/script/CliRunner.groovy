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
     * Instantiate the runner object creating a new session
     */
    def CliRunner( ) {
        this.session = new Session()
    }

    /**
     * Instantiate the runner object using the specified session
     *
     * @param session
     */
    def CliRunner(Session session) {
        this.session = session
    }

    /**
     * @return The interpreted script object
     */
    AbstractScript getScript() { script }

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
        if( !name ) {
            name = FilenameUtils.getBaseName(scriptFile.toString())
        }

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

        parseScript( file.text, args )
    }

    protected AbstractScript parseScript( String flow, String... args = null) {

        params['__$session'] = session
        params['args'] = args


        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addStaticStars(Nextflow.getName())

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = AbstractScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(TaskScriptClosureTransform))

        // run and wait for termination
        def groovy = new GroovyShell(this.class.classLoader, params, config)
        groovy.parse( flow ) as AbstractScript
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
        assert script, "Missing script to run"

        // -- define work-directory
        defineWorkDirectory()

        // -- launch the script execution
        def result = script.run()

        // -- normalize the
        if ( result == null ) {
            return null
        }
        else if( result instanceof Collection || result.class.isArray()) {
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
            println getVersion()
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


        File scriptFile = new File(options.arguments[0])
        String scriptName = FilenameUtils.getBaseName(options.arguments[0])
        if ( !scriptFile.exists() ) {
            System.err.println "The specified script file does not exist: '$scriptFile'"
            System.exit(1)
        }

        try {
            // -- create a new runner instance
            def runner = new CliRunner()

            // -- define the working directory
            runner.workDirectory = validateWorkDirectory( options.workDirectory, scriptName )
            log.info "Work directory: ${runner.workDirectory}"

            // -- define the variable bindings
            if ( options.params ) {
                options.params.each { k,v -> runner.params.setParam(k,v) }
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



    static String getVersion() {
        return Const.APP_VER
    }

}
