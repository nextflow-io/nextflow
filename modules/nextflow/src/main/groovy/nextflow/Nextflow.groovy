/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow

import static nextflow.file.FileHelper.*

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.ast.OpXform
import nextflow.ast.OpXformImpl
import nextflow.exception.StopSplitIterationException
import nextflow.exception.WorkflowScriptErrorException
import nextflow.extension.GroupKey
import nextflow.extension.OperatorImpl
import nextflow.file.FileHelper
import nextflow.file.FilePatternSplitter
import nextflow.mail.Mailer
import nextflow.script.TokenBranchDef
import nextflow.script.TokenMultiMapDef
import nextflow.splitter.FastaSplitter
import nextflow.splitter.FastqSplitter
import nextflow.util.ArrayTuple
import nextflow.util.CacheHelper
import org.slf4j.Logger
import org.slf4j.LoggerFactory
/**
 * Defines the main methods imported by default in the script scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Nextflow {

    static private Session getSession() { Global.session as Session }

    // note: groovy `Slf4j` annotation causes a bizarre issue
    // https://issues.apache.org/jira/browse/GROOVY-7371
    // declare public so that can be accessed from the user script
    public static final Logger log = LoggerFactory.getLogger(Nextflow)

    private static final Random random = new Random()


    static private fileNamePattern( FilePatternSplitter splitter, Map opts ) {

        final scheme = splitter.scheme
        final target = scheme ? "$scheme://$splitter.parent" : splitter.parent
        final folder = toCanonicalPath(target)
        final pattern = splitter.fileName

        if( opts == null ) opts = [:]
        if( !opts.type ) opts.type = 'file'

        def result = new LinkedList()
        try {
            FileHelper.visitFiles(opts, folder, pattern) { Path it -> result.add(it) }
        }
        catch (NoSuchFileException e) {
            log.debug "No such file or directory: $folder -- Skipping visit"
        }
        return result

    }

    static private String str0(value) {
        if( value==null )
            return null
        if( value instanceof CharSequence )
            return value.toString()
        if( value instanceof File )
            return value.toString()
        if( value instanceof Path )
            return value.toUriString()
        throw new IllegalArgumentException("Invalid file path type - offending value: $value [${value.getClass().getName()}]")
    }

    /**
     * Get one or more file object given the specified path or glob pattern.
     *
     * @param options
     *      A map holding one or more of the following options:
     *      - type: Type of paths returned, either `file`, `dir` or `any` (default: `file`)
     *      - hidden: When `true` includes hidden files in the resulting paths (default: `false`)
     *      - maxDepth: Maximum number of directory levels to visit (default: no limit)
     *      - followLinks: When `true` it follows symbolic links during directories tree traversal, otherwise they are managed as files (default: `true`)
     *
     * @param path A file path eventually including a glob pattern e.g. /some/path/file*.txt
     * @return An instance of {@link Path} when a single file is matched or a list of {@link Path}s
     */
    static file( Map options = null, def filePattern ) {

        if( !filePattern )
            throw new IllegalArgumentException("Argument of `file` function cannot be ${filePattern==null?'null':'empty'}")

        final path = filePattern as Path
        final glob = options?.containsKey('glob') ? options.glob as boolean : isGlobAllowed(path)
        if( !glob ) {
            return checkIfExists(path, options)
        }

        // if it isn't a glob pattern simply return it a normalized absolute Path object
        final strPattern = str0(filePattern)
        final splitter = FilePatternSplitter.glob().parse(strPattern)
        if( !splitter.isPattern() ) {
            final normalised = splitter.strip(strPattern)
            return checkIfExists(asPath(normalised), options)
        }

        // revolve the glob pattern returning all matches
        return fileNamePattern(splitter, options)
    }

    static files( Map options=null, def path ) {
        def result = file(options, path)
        return result instanceof List ? result : [result]
    }


    /**
     * Creates a {@link ArrayTuple} object with the given open array items
     *
     * @param args The items used to created the tuple object
     * @return An instance of {@link ArrayTuple} populated with the given argument(s)
     */
    static  ArrayTuple tuple( def value ) {
        if( !value )
            return new ArrayTuple()

        new ArrayTuple( value instanceof Collection ? (Collection)value : [value] )
    }

    /**
     * Creates a {@link ArrayTuple} object with the given open array items
     *
     * @param args The items used to created the tuple object
     * @return An instance of {@link ArrayTuple} populated with the given argument(s)
     */
    static ArrayTuple tuple( Object ... args ) {
        tuple( args as List )
    }

    /**
     * Creates a FASTQ splitter handler for the given object
     *
     * @param obj The object to be managed as a FASTQ
     * @return An instance of {@link FastqSplitter
     */
    @Deprecated
    static FastqSplitter fastq( obj ) {
        (FastqSplitter)new FastqSplitter('fastq').target(obj)
    }

    /**
     * Creates a FASTA splitter handler for the given object
     *
     * @param obj The object to be managed as a FASTA
     * @return An instance of {@link FastqSplitter
     */
    @Deprecated
    static FastaSplitter fasta( obj ) {
        (FastaSplitter)new FastaSplitter('fasta').target(obj)
    }

    /**
     * Interrupts the iteration when using a split operators
     */
    static void stop() {
        throw new StopSplitIterationException()
    }


    /**
     * Stop the current execution returning an error code and message
     *
     * @param exitCode The exit code to be returned
     * @param message The message that will be reported in the log file (optional)
     */
    @Deprecated
    static void exit(int exitCode, String message = null) {
        if( session.aborted ) {
            log.debug "Ignoring exit because execution is already aborted -- message=$message"
            return
        }
        
        if ( exitCode && message ) {
            log.error message
        }
        else if ( message ) {
            log.info message
        }
        System.exit(exitCode)
    }

    /**
     * Stop the current execution returning a 0 error code and the specified message
     *
     * @param message The message that will be reported in the log file
     */
    @Deprecated
    static void exit( String message ) {
        exit(0, message)
    }

    /**
     * Throws a script runtime error
     * @param message An optional error message
     */
    static void error( String message = null ) {
        throw message ? new WorkflowScriptErrorException(message) : new WorkflowScriptErrorException()
    }

    /**
     * Create a folder for the given key. It guarantees to return the same folder name
     * the same provided object key.
     *
     * @param key An object to be used as cache-key creating the folder, it can be any object
     *          or an array or objects to use multi-objects key
     *
     * @return The {@code Path} to the cached directory or a newly created folder for the specified key
     */
    @Deprecated
    static Path cacheableDir( Object key ) {
        assert key, "Please specify the 'key' argument on 'cacheableDir' method"

        final session = (Session)Global.session
        if( !session )
            throw new IllegalStateException("Invalid access to `cacheableDir` method -- Session object not yet defined")

        def hash = CacheHelper.hasher([ session.uniqueId, key, session.cacheable ? 0 : random.nextInt() ]).hash()

        def file = FileHelper.getWorkFolder(session.workDir, hash)
        if( !file.exists() && !file.mkdirs() ) {
            throw new IOException("Unable to create directory: $file -- Check file system permissions" )
        }

        return file
    }

    /**
     * Create a file for the given key. It guarantees to return the same file name
     * the same provided object key.
     *
     * @param key
     * @param name
     * @return
     */
    @Deprecated
    static Path cacheableFile( Object key, String name = null ) {

        // the cacheability is guaranteed by the folder
        def folder = cacheableDir(key)

        if( !name ) {
            if( key instanceof File ) {
                name =  key.getName()
            }
            else if( key instanceof Path ) {
                name =  key.getName()
            }
            else {
                name = key.toString()
            }
        }

        return folder.resolve(name)
    }

    /**
     * This method is exposed as a public API to script, it should be removed
     *
     * @return Create a temporary directory
     */
    @Deprecated
    static Path tempDir( String name = null, boolean create = true ) {
        final session = (ISession)Global.session
        if( !session )
            throw new IllegalStateException("Invalid access to `tempDir` method -- Session object not yet defined")

        def path = FileHelper.createTempFolder(session.workDir)
        if( name )
            path = path.resolve(name)

        if( !path.exists() && create && !path.mkdirs() )
            throw new IOException("Unable to create directory: $path -- Check file system permissions" )

        return path
    }

    /**
     * This method is exposed as a public API to script, it should be removed
     *
     * @return Create a temporary file
     */
    @Deprecated
    static Path tempFile( String name = null, boolean create = false ) {

        if( !name )
            name = 'file.tmp'

        def folder = tempDir()
        def result = folder.resolve(name)
        if( create )
            Files.createFile(result)

        return result
    }

    /**
     * Implements built-in send mail functionality
     *
     * @param params
     *      A map object holding the mail parameters
     *      - to: String or a List of strings representing the mail recipients
     *      - cc: String or a List of strings representing the mail recipients in carbon copy
     *      - from: String representing the mail sender address
     *      - subject: The mail subject
     *      - content: The mail content
     *      - attach: One or more list attachment
     */
    static void sendMail( Map params ) {

        new Mailer()
            .setConfig(Global.session.config.mail as Map)
            .send(params)

    }

    /**
     * Implements built-in send mail functionality
     *
     * @param params
     *    A closure representing the mail message to send eg
     *    <code>
     *        sendMail {
     *          to 'me@dot.com'
     *          from 'your@name.com'
     *          attach '/some/file/path'
     *          subject 'Hello'
     *          content '''
     *           Hi there,
     *           Hope this email find you well
     *          '''
     *        }
     *    <code>
     */
    static void sendMail( Closure params ) {
        new Mailer()
                .setConfig(Global.session.config.mail as Map)
                .send(params)
    }

    /**
     * Creates a groupTuple dynamic key
     *
     * @param key
     * @param size
     * @return
     */
    static GroupKey groupKey(key, int size) {
        new GroupKey(key,size)
    }

    /**
     * Marker method to create a closure to be passed to {@link OperatorImpl#branch(DataflowReadChannel, groovy.lang.Closure)}
     * operator.
     *
     * Despite apparently is doing nothing, this method is needed as marker to apply the {@link OpXform} AST
     * transformation required to interpret the closure content as required for the branch evaluation.
     *
     * @see OperatorImpl#branch(DataflowReadChannel, Closure)
     * @see OpXformImpl
     *
     * @param closure
     * @return
     */
    static Closure<TokenBranchDef> branchCriteria(Closure<TokenBranchDef> closure) { closure }

    /**
     * Marker method to create a closure to be passed to {@link OperatorImpl#fork(DataflowReadChannel, Closure)}
     * operator.
     *
     * Despite apparently is doing nothing, this method is needed as marker to apply the {@link OpXform} AST
     * transformation required to interpret the closure content as required for the branch evaluation.
     *
     * @see OperatorImpl#multiMap(groovyx.gpars.dataflow.DataflowReadChannel, groovy.lang.Closure) (DataflowReadChannel, Closure)
     * @see OpXformImpl
     *
     * @param closure
     * @return
     */
    static Closure<TokenMultiMapDef> multiMapCriteria(Closure<TokenBranchDef> closure) { closure }

}
