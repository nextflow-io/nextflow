/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.script

import java.nio.file.Path

import com.google.common.hash.Hashing
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Channel
import nextflow.Const
import nextflow.Nextflow
import nextflow.Session
import nextflow.ast.NextflowDSL
import nextflow.ast.NextflowXform
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.file.FileHelper
import nextflow.util.ConfigHelper
import nextflow.util.HistoryFile
import nextflow.util.VersionNumber
import org.apache.commons.lang.StringUtils
import org.codehaus.groovy.control.CompilerConfiguration
import org.codehaus.groovy.control.customizers.ASTTransformationCustomizer
import org.codehaus.groovy.control.customizers.ImportCustomizer
import static nextflow.util.ConfigHelper.parseValue
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
     * The pipeline file (it may be null when it's provided as string)
     */
    private ScriptFile scriptFile

    /**
     * The pipeline as a text content
     */
    private String scriptText

    /*
     * the script raw output
     */
    private def output

    /**
     * The script result
     */
    private def result

    /**
     * The launcher command line
     */
    private String commandLine

    /**
     * The used configuration profile
     */
    String profile

    /**
     * Instantiate the runner object creating a new session
     */
    ScriptRunner( ) {
        this( [:] )
    }

    ScriptRunner( Map config ) {
        this.session = new Session(config)
    }

    ScriptRunner( Session session ) {
        this.session = session
    }

    ScriptRunner setScript( ScriptFile script ) {
        this.scriptFile = script
        this.scriptText = script.text
        return this
    }

    ScriptRunner( ConfigBuilder builder ) {
        this.session = new Session(builder.build())
        // note config files are collected during the build process
        // this line should be after `ConfigBuilder#build`
        this.session.configFiles = builder.parsedConfigFiles
    }

    ScriptRunner setScript( String text ) {
        this.scriptText = text
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
        assert scriptText

        // init session
        session.init(scriptFile?.main)

        // start session
        session.start()
        try {
            // parse the script
            script = parseScript(scriptText, args)
            // validate the config
            validate()
            // run the code
            run()
            // await termination
            terminate()
        }
        catch (Throwable e) {
            session.abort(e)
            throw e
        }

        if( !session.success ) {
            throw new AbortRunException()
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
        assert scriptText
        assert methodName

        // init session
        session.init(scriptFile.main)

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
        session.binding.setArgs( new ArgsList(args) )
        session.binding.setParams( (Map)session.config.params )
        // TODO add test for this property
        session.binding.setVariable( 'baseDir', session.baseDir )
        session.binding.setVariable( 'workDir', session.workDir )
        if( scriptFile ) {
            def meta = new WorkflowMetadata(this)
            session.binding.setVariable( 'workflow', meta )
            session.binding.setVariable( 'nextflow', meta.nextflow )
        }

        // generate an unique class name
        session.scriptClassName = generateClassName(scriptText)

        // define the imports
        def importCustomizer = new ImportCustomizer()
        importCustomizer.addImports( StringUtils.name, groovy.transform.Field.name )
        importCustomizer.addImports( Path.name )
        importCustomizer.addImports( Channel.name )
        importCustomizer.addStaticStars( Nextflow.name )

        def config = new CompilerConfiguration()
        config.addCompilationCustomizers( importCustomizer )
        config.scriptBaseClass = BaseScript.class.name
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowDSL))
        config.addCompilationCustomizers( new ASTTransformationCustomizer(NextflowXform))

        // extend the class-loader if required
        def gcl = new GroovyClassLoader()
        def libraries = ConfigHelper.resolveClassPaths( session.getLibDir() )

        libraries?.each { Path lib -> def path = lib.complete()
            log.debug "Adding to the classpath library: ${path}"
            gcl.addClasspath(path.toString())
        }

        // set the byte-code target directory
        def targetDir = FileHelper.createLocalDir()
        config.setTargetDirectory(targetDir.toFile())
        // add the directory of generated classes to the lib path
        // so that it can be propagated to remote note (when necessary)
        session.getLibDir().add(targetDir)

        // set the script class-loader
        session.classLoader = gcl

        // run and wait for termination
        BaseScript result
        def groovy = new GroovyShell(gcl, session.binding, config)
        result = groovy.parse(scriptText, session.scriptClassName) as BaseScript

        session.onShutdown { targetDir.deleteDir() }
        return result
    }

    /**
     * Creates a unique name for the main script class in order to avoid collision
     * with the implicit and user variables
     */
    @PackageScope String generateClassName(String text) {
        def hash = Hashing
                .murmur3_32()
                .newHasher()
                .putUnencodedChars(text)
                .hash()

        return "_nf_script_${hash}"
    }

    /**
     * Check preconditions before run the main script
     */
    protected void validate() {
        checkConfig()
        checkVersion()
    }

    @PackageScope void checkConfig() {
        session.validateConfig(script.getProcessNames())
    }

    @PackageScope VersionNumber getCurrentVersion() {
        new VersionNumber(Const.APP_VER)
    }

    @PackageScope void checkVersion() {
        def version = session.manifest.getNextflowVersion()?.trim()
        if( !version )
            return

        // when the version string is prefix with a `!`
        // an exception is thrown is the version does not match
        boolean important = false
        if( version.startsWith('!') ) {
            important = true
            version = version.substring(1).trim()
        }

        if( !getCurrentVersion().matches(version) ) {
            important ? showVersionError(version) : showVersionWarning(version)
        }
    }

    @PackageScope void showVersionError(String ver) {
        throw new AbortOperationException("Nextflow version $Const.APP_VER does not match workflow required version: $ver")
    }

    @PackageScope void showVersionWarning(String ver) {
        log.warn "Nextflow version $Const.APP_VER does not match workflow required version: $ver -- Execution will continue, but things may break!"
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
        session.cleanup()
        log.debug "> Execution complete -- Goodbye"
    }

    /**
     * @param cli The launcher command line string
     */
    void verifyAndTrackHistory(String cli, String name) {

        // -- when resume, make sure the session id exists in the executions history
        if( session.resumeMode && !HistoryFile.DEFAULT.checkExistsById(session.uniqueId.toString())) {
            throw new AbortOperationException("Can't find a run with the specified id: ${session.uniqueId} -- Execution can't be resumed")
        }

        if( !cli )
            return
        def p = cli.indexOf('nextflow ')
        commandLine = p != -1 ? 'nextflow ' + cli.substring(p+9) : cli
        def revisionId = scriptFile.commitId ?: scriptFile.scriptId
        HistoryFile.DEFAULT.write( name, session.uniqueId, revisionId, commandLine )
    }


    @PackageScope
    ScriptFile getScriptFile() { scriptFile }

    @PackageScope
    String getCommandLine() { commandLine }

    @PackageScope
    @CompileDynamic
    def fetchContainers() {

        def result = [:]
        if( session.config.process instanceof Map<String,?> ) {

            /*
             * look for `container` definition at process level
             */
            session.config.process.each { String name, value ->
                if( name.startsWith('$') && value instanceof Map && value.container ) {
                    result[name] = resolveClosure(value.container)
                }
            }

            /*
             * default container definition
             */
            def container = session.config.process.container
            if( container ) {
                if( result ) {
                    result['default'] = resolveClosure(container)
                }
                else {
                    result = resolveClosure(container)
                }
            }

        }

        return result
    }

    /**
     * Resolve dynamically defined attributes to the actual value
     *
     * @param val A process container definition either a plain string or a closure
     * @return The actual container value
     */
    protected String resolveClosure( val ) {
        if( val instanceof Closure ) {
            try {
                return val.cloneWith(session.binding).call()
            }
            catch( Exception e ) {
                log.debug "Unable to resolve dynamic `container` directive -- cause: ${e.message ?: e}"
                return "(dynamic resolved)"
            }
        }

        return String.valueOf(val)
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
