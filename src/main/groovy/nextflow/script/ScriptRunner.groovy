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

import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.cli.CmdRun
import nextflow.exception.AbortOperationException
import nextflow.util.ConfigHelper
import org.apache.commons.lang.StringUtils
import org.apache.commons.lang.exception.ExceptionUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
/**
 * Application main class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ScriptRunner {

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
    private ScriptBinding bindings

    /**
     * The pipeline file (it may be null when it's provided as string)
     */
    private Path scriptFile

    /**
     * The pipeline as a text content
     */
    private String setScript

    /**
     * The command line options
     */
    private CmdRun runOptions

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
    def ScriptRunner( ) {
        this( [:] )
    }

    def ScriptRunner( Map config ) {
        session = new Session(config)
        bindings = new ScriptBinding(session)
    }

    def ScriptRunner setScript( Path file ) {
        this.scriptFile = file
        this.setScript = file.text
        return this
    }

    def ScriptRunner setScript( String text ) {
        this.setScript = text
        return this
    }

    def ScriptRunner setRunOptions( CmdRun options ) {
        this.runOptions = options
        return this
    }

    Session getSession() { session }

    /**
     * @return The interpreted script object
     */
    BaseScript getScriptObj() { script }

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

    def execute( List<String> args = null ) {
        assert setScript

        // init session
        if( runOptions )
            session.init(runOptions, scriptFile)

        // start session
        session.start()
        try {
            // parse the script
            script = parseScript(setScript, args)
            // run the code
            run()
        }
        catch( MissingPropertyException e ) {
            throw new AbortOperationException(getErrorMessage(e, session.scriptName), e)
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
    def test ( String methodName, List<String> args = null ) {
        assert setScript
        assert methodName

        // init session
        session.init(runOptions, scriptFile)

        script = parseScript(setScript, args)
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


    def normalizeOutput() {
        if( output instanceof Object[] ) {
            result = output as Collection
        }
        else {
            result = output
        }
    }

    protected BaseScript parseScript( File file, List<String> args = null ) {
        assert file
        parseScript( file.text, args )
    }

    protected BaseScript parseScript( String scriptText, List<String> args = null) {
        log.debug "> Script parsing"
        bindings.setArgs( new ArgsList(args) )
        bindings.setParams( session.config.params as Map )
        // TODO add test for this property
        bindings.setVariable( 'baseDir', session.baseDir )
        bindings.setVariable( 'workDir', session.workDir )

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
        def libraries = ConfigHelper.resolveClassPaths( session.getLibDir() )

        libraries?.each { Path lib -> def path = lib.fixed()
            log.debug "Adding to the classpath library: ${path}"
            gcl.addClasspath(path.toString())
        }

        // set the byte-code target directory
        def targetDir = Files.createTempDirectory('nxf')
        config.setTargetDirectory(targetDir.toFile())
        // add the directory of generated classes to the lib path
        // so that it can be propagated to remote note (when necessary)
        session.getLibDir().add(targetDir)

        // run and wait for termination
        BaseScript result
        def groovy = new GroovyShell(gcl, bindings, config)
        if( session.scriptName )
            result = groovy.parse(scriptText, session.scriptName) as BaseScript
        else
            result = groovy.parse(scriptText) as BaseScript

        session.onShutdown { targetDir.deleteDir() }
        return result
    }


    /**
     * Launch the Nextflow script execution
     *
     * @return The value as returned by the user provided script
     */
    protected run() {
        log.debug "> Launching execution"
        assert script, "Missing script instance to run"

        // -- launch the script execution
        output = script.run()
    }

    protected terminate() {
        log.debug "> Await termination "
        session.await()
        normalizeOutput()
        session.destroy()
    }


    /**
     * Find out the script line where the error has thrown
     */
    static String getErrorMessage( Throwable e, String scriptName ) {

        def lines = ExceptionUtils.getStackTrace(e).split('\n')
        List error = null
        for( String str : lines ) {
            if( (error=getErrorLine(str,scriptName))) {
                break
            }
        }

        if( e instanceof MissingPropertyException ) {
            def result = "Not such variable: '${e.getProperty()}'"
            if( error )
                result += " -- check script '${error[0]}' at line: ${error[1]}"
            else if( scriptName )
                result += " -- check that you don't mispelled it in scipt '$scriptName'"

            return result
        }
        else if( error ) {
            (e.message ?: e.toString()) + " at $error"
        }
        else {
            e.message ?: e.toString()
        }
    }


    @PackageScope
    static List<String> getErrorLine( String line, String scriptName = null) {
        if( scriptName==null )
            scriptName = '.+'

        def pattern = ~/.*\(($scriptName):(\d*)\).*/
        def m = pattern.matcher(line)
        return m.matches() ? [m.group(1), m.group(2)] : null
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
                throw new AbortOperationException("Argument array index cannot be lower than zero")
            }

            if( pos >= size() ) {
                throw new AbortOperationException("Arguments index out of range: $pos -- You may have not entered all arguments required by the pipeline")
            }

            super.get(pos)
        }

        def String getAt(int pos) {
            get(pos)
        }

    }

}
