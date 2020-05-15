/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.config.ConfigBuilder
import nextflow.exception.AbortOperationException
import nextflow.exception.AbortRunException
import nextflow.util.HistoryFile
import static nextflow.util.ConfigHelper.parseValue
/**
 * Run a nextflow script file
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
     * The script interpreter
     */
    private ScriptParser scriptParser

    /**
     * The pipeline file (it may be null when it's provided as string)
     */
    private ScriptFile scriptFile

    /**
     * The script result
     */
    private def result

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

    ScriptRunner setScript( Path script ) {
        setScript(new ScriptFile(script))
    }

    ScriptRunner setScript( ScriptFile script ) {
        this.scriptFile = script
        return this
    }

    ScriptRunner( ConfigBuilder builder ) {
        this.session = new Session(builder.build())
        // note config files are collected during the build process
        // this line should be after `ConfigBuilder#build`
        this.session.configFiles = builder.parsedConfigFiles
    }

    Session getSession() { session }

    /**
     * @return The interpreted script object
     */
    @Deprecated BaseScript getScriptObj() { scriptParser.script }

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

    def execute( List<String> args = null, String entryName=null ) {
        assert scriptFile

        // init session
        session.init(scriptFile, args)

        // start session
        session.start()
        try {
            // parse the script
            parseScript(scriptFile, entryName)
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
        assert methodName
        // init session
        session.init(scriptFile, args)
        parseScript(scriptFile, null)
        def values = args ? args.collect { parseValue(it) } : null

        def methodsToTest
        if ( methodName == '%all' ) {
            methodsToTest = scriptObj.metaClass.getMethods().findAll { it.name.startsWith('test') }.collect { it.name }.unique()
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
                def metaMethod = scriptObj.metaClass.getMethods().find { it.name == name }

                if( !metaMethod ) {
                    log.error "Unknown function '$name' -- Check the name spelling"
                }
                else {
                    def result = metaMethod.invoke(scriptObj, values as Object[] )
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

    def normalizeOutput(output) {
        return output
    }

    protected void parseScript( ScriptFile scriptFile, String entryName ) {
        scriptParser = new ScriptParser(session)
                            .setEntryName(entryName)
                            .parse(scriptFile.main)
        session.script = scriptParser.script
    }


    /**
     * Launch the Nextflow script execution
     *
     * @return The value as returned by the user provided script
     */
    protected run() {
        log.debug "> Launching execution"
        assert scriptParser, "Missing script instance to run"
        // -- launch the script execution
        scriptParser.runScript()
        // -- normalise output
        result = normalizeOutput(scriptParser.getResult())
        // -- ignite dataflow network
        session.fireDataflowNetwork()
    }

    protected terminate() {
        log.debug "> Await termination "
        session.await()
        session.destroy()
        session.cleanup()
        log.debug "> Execution complete -- Goodbye"
    }

    /**
     * @param cli The launcher command line string
     */
    void verifyAndTrackHistory(String cli, String name) {
        assert cli, 'Missing launch command line'

        final ignore = System.getenv('NXF_IGNORE_RESUME_HISTORY') as Boolean
        // -- when resume, make sure the session id exists in the executions history
        if( session.resumeMode && !ignore && !HistoryFile.DEFAULT.checkExistsById(session.uniqueId.toString()) ) {
            throw new AbortOperationException("Can't find a run with the specified id: ${session.uniqueId} -- Execution can't be resumed")
        }

        def revisionId = scriptFile.commitId ?: scriptFile.scriptId
        HistoryFile.DEFAULT.write( name, session.uniqueId, revisionId, cli )
    }


    @PackageScope
    ScriptFile getScriptFile() { scriptFile }

    /**
     * Extends an {@code ArrayList} class adding a nicer index-out-of-range error message
     */
    static class ArgsList extends ArrayList<String> {


        ArgsList(List<String> values) {
            super( values ?: [] )
        }

        String get( int pos ) {
            if( pos < 0 ) {
                throw new AbortOperationException("Argument array index cannot be lower than zero")
            }

            if( pos >= size() ) {
                throw new AbortOperationException("Arguments index out of range: $pos -- You may have not entered all arguments required by the pipeline")
            }

            super.get(pos)
        }

        String getAt(int pos) {
            get(pos)
        }

    }

}
