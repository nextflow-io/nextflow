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

import java.nio.file.Path

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.exception.AbortOperationException
import nextflow.exception.MissingLibraryException
import nextflow.util.FileHelper
import org.apache.commons.io.FilenameUtils
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
     * The script file
     */
    private File scriptFile

    private String scriptName

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
    def ScriptRunner( ) {
        this( [:] )
    }

    def ScriptRunner( Map config ) {
        session = new Session(config)
        bindings = new ScriptBinding(session)
    }

    /**
     * Only for testing purposes
     *
     * @param config
     */
    def ScriptRunner( ConfigObject config ) {
        this(ConfigBuilder.configToMap(config))
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
        this.scriptName = FilenameUtils.getBaseName(scriptFile.toString())

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
        catch( MissingPropertyException e ) {
            throw new RuntimeException(getErrorMessage(e, scriptName), e)
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
        this.scriptName = FilenameUtils.getBaseName(scriptFile.toString())
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
        // TODO add test for this property
        bindings.setVariable( 'scriptDir', session.baseDir?.toPath() )

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

        // set the byte-code target directory
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
     * Find out the script line where the error has thrown
     */
    static getErrorMessage( Throwable e, String scriptName ) {

        def lines = ExceptionUtils.getStackTrace(e).split('\n')
        def error = null
        for( String str : lines ) {
            if( (error=getErrorLine(str,scriptName))) {
                break
            }
        }

        return (e.message ?: e.toString()) + ( error ? " at $error" : '' )
    }


    @PackageScope
    static String getErrorLine( String line, String scriptName = null) {
        if( scriptName==null )
            scriptName = '.+'

        def pattern = ~/.*\(($scriptName\.nf:\d*)\).*/
        def m = pattern.matcher(line)
        return m.matches() ? m.group(1) : null
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
