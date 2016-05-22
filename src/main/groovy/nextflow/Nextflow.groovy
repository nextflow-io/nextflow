/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

package nextflow

import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.Path

import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.exception.ProcessNotRecoverableException
import nextflow.exception.StopSplitIterationException
import nextflow.file.FileHelper
import nextflow.splitter.FastaSplitter
import nextflow.splitter.FastqSplitter
import nextflow.util.ArrayTuple
import nextflow.util.CacheHelper
/**
 * Defines the main methods imported by default in the script scope
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class Nextflow {

    static private final Random random = new Random()

    /**
     * Create a {@code DataflowVariable} binding it to the specified value
     *
     * @param value
     * @return
     */
    static <T> DataflowVariable<T> variable( T value = null ) {
        def result = new DataflowVariable<T>()
        if( value != null ) {
            result.bind(value)
        }
        result
    }

    /**
     * Create a {@code DataflowQueue} populating with the specified values
     * <p>
     * This 'queue' data structure can be viewed as a point-to-point (1 to 1, many to 1) communication channel.
     * It allows one or more producers send messages to one reader.
     *
     * @param values
     * @return
     */
    static DataflowQueue channel( Collection values = null ) {

        final channel = new DataflowQueue()
        if ( values )  {
            // bind e
            values.each { channel.bind(it)  }

            // since queue is 'finite' close it by a poison pill
            // so the operator will stop on when all values in the queue are consumed
            // (otherwise it will wait forever for a new entry)
            channel << Channel.STOP
        }

        return channel
    }

    /**
     * Create a {@code DataflowQueue} populating with a single value
     * <p>
     * This 'queue' data structure can be viewed as a point-to-point (1 to 1, many to 1) communication channel.
     * It allows one or more producers send messages to one reader.
     *
     * @param item
     * @return
     */
    static DataflowQueue channel( Object... items ) {
        return channel(items as List)
    }

    static private fileNamePattern( String path, Map opts, java.nio.file.FileSystem fs ) {

        def ( String folder, String pattern, String scheme ) = FileHelper.getFolderAndPattern(path)
        if( !fs )
            fs = FileHelper.fileSystemForScheme(scheme)

        if( opts == null ) opts = [:]
        if( !opts.type ) opts.type = 'file'

        def result = new LinkedList()
        try {
            FileHelper.visitFiles(opts, fs.getPath(folder), pattern) { Path it -> result.add(it) }
        }
        catch (NoSuchFileException e) {
            log.debug "No such file: $folder -- Skipping visit"
        }
        return result

    }

    static file( Map options = null, def path ) {
        assert path

        if( path == null )
            return null

        final isPath = path instanceof Path
        final str = path.toString()
        final fs = isPath ? path.getFileSystem() : null

        // if it isn't a glob pattern simply return it a normalized absolute Path object
        if( !FileHelper.isGlobPattern(str) ) {
            if( !isPath ) path = FileHelper.asPath(str)
            return path.complete()
        }

        // revolve the glob pattern returning all matches
        return fileNamePattern(path.toString(), options, fs)
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
    static FastqSplitter fastq( obj ) {
        (FastqSplitter)new FastqSplitter('fastq').target(obj)
    }

    /**
     * Creates a FASTA splitter handler for the given object
     *
     * @param obj The object to be managed as a FASTA
     * @return An instance of {@link FastqSplitter
     */
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
    static void exit(int exitCode, String message = null) {
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
    static void exit( String message ) {
        exit(0, message)
    }

    /**
     * Throws a script runtime error
     * @param message An optional error message
     */
    static void error( String message = null ) {
        throw message ? new ProcessNotRecoverableException(message) : new ProcessNotRecoverableException()
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
    static Path cacheableDir( Object key ) {
        assert key, "Please specify the 'key' argument on 'cacheableDir' method"

        final session = (Session)Global.session
        if( !session )
            throw new IllegalStateException("Invalid access to `cacheableDir` method -- Session object not yet defined")

        def hash = CacheHelper.hasher([ session.uniqueId, key, session.cacheable ? 0 : random.nextInt() ]).hash()

        def file = FileHelper.getWorkFolder(session.workDir, hash)
        if( !file.exists() && !file.mkdirs() ) {
            throw new IOException("Unable to create folder: $file -- Check file system permission" )
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
     * @return Create a temporary directory
     */
    static Path tempDir( String name = null, boolean create = true ) {
        final session = (ISession)Global.session
        if( !session )
            throw new IllegalStateException("Invalid access to `tempDir` method -- Session object not yet defined")

        def path = FileHelper.createTempFolder(session.workDir)
        if( name )
            path = path.resolve(name)

        if( !path.exists() && create && !path.mkdirs() )
            throw new IOException("Unable to create folder: $path -- Check file system permission" )

        return path
    }

    /**
     * @return Create a temporary file
     */
    static Path tempFile( String name = null, boolean create = false ) {

        if( !name )
            name = 'file.tmp'

        def folder = tempDir()
        def result = folder.resolve(name)
        if( create )
            Files.createFile(result)

        return result
    }



}
