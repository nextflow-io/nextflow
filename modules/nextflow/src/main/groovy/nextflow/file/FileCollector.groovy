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

package nextflow.file
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.NoSuchFileException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx
import nextflow.io.SkipLinesInputStream
import nextflow.util.CacheHelper
import nextflow.util.KryoHelper
/**
 * File collector base class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class FileCollector implements Closeable {

    static final public OpenOption[] APPEND = [StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.APPEND] as OpenOption[]

    static final public OpenOption[] TRUNCATE = [StandardOpenOption.CREATE, StandardOpenOption.WRITE, StandardOpenOption.TRUNCATE_EXISTING] as OpenOption[]

    Boolean newLine

    def seed

    Integer skipLines

    boolean keepHeader

    boolean deleteTempFilesOnClose = true;

    CacheHelper.HashMode hashMode = CacheHelper.HashMode.STANDARD

    List hashKeys

    boolean cacheable

    boolean resumable

    private Path tempDir

    private boolean created

    private boolean closed

    protected FileCollector() { }

    /**
     * Create the temporary directory
     * @return A {@code Path} instance to the directory location
     */
    protected Path createTempDir() {

        if( !tempDir ) {
            return tempDir = FileHelper.createLocalDir()
        }

        // read the path attributes
        def attr = Files.readAttributes(tempDir, BasicFileAttributes)
        // when it is a dir => create a temp file there
        if ( attr.isDirectory() ) {
            return tempDir
        }

        // if its a file use it
        if( attr.isRegularFile()) {
            throw new FileAlreadyExistsException("Cannot create sort path: $tempDir -- A file with the same name already exists")
        }

        return Files.createDirectories(tempDir)
    }

    /**
     * @param path the temporary directory to be used
     */
    public void setTempDir( Path path ) {
        this.tempDir = path
        created = false
    }

    /**
     * Retrieve the collector temporary folder i.e. working directory. The first time, it invokes {@link #createTempDir()} method
     * to create the directory path
     *
     * @return Collector temporary directory
     */
    Path getTempDir() {
        if( !created ) {
            tempDir = createTempDir()
            created = true
        }
        return tempDir
    }

    /**
     * Add and entry to the collector set
     *
     * @param key A grouping key. All entries with the same key will belong in the same grouping set.
     * @param value A value to be added
     * @return The collector {@code FileCollector} object itself
     */
    abstract FileCollector add( String key, value );

    /**
     * Save an collect entry to the respective grouping file
     * @param closure Closure mapping a group name to the respective grouping file e.g. { name -> Paths.get(name) }
     */
    abstract void saveFile( Closure<Path> closure );

    protected HashCode makeHash() { return null }

    /**
     * Normalize values to a {@link InputStream}
     *
     * @param value The user provided value
     * @return An {@link InputStream} referring the value
     */
    protected InputStream normalizeToStream0( value ) {
        if( value instanceof Path )
            return value.newInputStream()

        if( value instanceof File )
            return value.newInputStream()

        if( value instanceof CharSequence )
            return new ByteArrayInputStream(value.toString().getBytes())

        if( value instanceof byte[] )
            return new ByteArrayInputStream(value as byte[])

        throw new IllegalArgumentException("Not a valid file collector argument [${value.class.name}]: $value")
    }

    protected InputStream normalizeToStream( value ) {
        def result = normalizeToStream0(value)
        if( result && skipLines ) {
            result = new SkipLinesInputStream(result,skipLines)
            if( keepHeader )
                result.consumeHeader()
        }
        return result
    }

    protected void appendHeader(InputStream data, Object name, OutputStream target) {
        assert !(keepHeader && seed), "Argument `keepHeader` and `seed` conflicts"

        def header
        if( keepHeader && data instanceof SkipLinesInputStream) {
            header = data.getHeader()
        }
        else {
            header = seed instanceof Closure ? ((Closure)((Closure)seed).clone()).call(name) : seed
        }

        InputStream result = null
        if( header instanceof Map && header.get(name) ) {
            result = normalizeToStream0(header.get(name))
        }
        else if ( header )
            result = normalizeToStream0(header)

        if( result )
            appendStream(result, target)
    }


    /**
     * Append the content of a file to the target file
     *
     * @param source The source stream representing the data to append to {@code  target}
     * @param target The target object
     */
    protected void appendStream( InputStream source, OutputStream target ) {
        int n
        byte[] buffer = new byte[10 * 1024]

        try {
            while( (n=source.read(buffer)) > 0 ) {
                target.write(buffer,0,n)
            }
            // append the new line separator
            if( newLine )
                target.write( System.lineSeparator().bytes )
        }
        finally {
            source.closeQuietly()
        }
    }


    /**
     * Save the entries collected grouping them into files whose name is given by the
     * correspondent group key
     *
     * @param target A {@link Path} to the folder where files will be saved. If the folder does not exists, is created
     *              automatically
     * @return The list of files where entries have been saved.
     */
    private List<Path> saveTo0(Path target) {
        target.createDirIfNotExists()

        def result = []
        saveFile { String name ->
            Path newFile = target.resolve(name)
            result << newFile
            return newFile
        }
        return result
    }

    List<Path> saveTo(Path target) {

        // verify if a cached list exists
        final hash = makeHash()?.toString()
        Path temp = hash ? FileHelper.getLocalTempPath().resolve("${hash}.collect-file") : null

        // try to retrieve cached files
        List<Path> items = null
        if( resumable && temp ) {
            items = retrieveCachedFiles(temp)
        }

        // get the list of files to generated
        if( items == null ) {
            items = saveTo0(target)
            // save the list of collected files
            if( cacheable && temp ) {
                cacheCollectedFile(items, temp)
            }
        }

        return items
    }

    /**
     * Close the collector deleting temporary files
     *
     * See {@link #deleteTempFilesOnClose} {@link #safeClose()}
     */
    void close() {
        if( tempDir && deleteTempFilesOnClose ) {
            log.debug "Deleting file collector temp dir: ${tempDir}"
            FilesEx.deleteDir(tempDir);
        }
    }

    /**
     * A synchronised version of {@code #close} operation
     */
    synchronized void safeClose() {
        if( closed ) return
        try {
            close()
            closed = true
        }
        catch( Exception e ) {
            log.debug("Unable to close file collector",e)
        }
    }

    private List<Path> retrieveCachedFiles(Path temp) {
        try {
            log.debug ">> temp file exists? ${temp.exists()}"
            List<Path> items = (List<Path>)KryoHelper.deserialize(temp)
            def notFound = items.find { Path p -> !p.exists() }
            if( notFound ) throw new NoSuchFileException("Missing cached file: $notFound")
            log.debug "Retrivied cached collect-files from: ${temp} -- cached files: ${items}"
            return items
        }
        catch( NoSuchFileException e ) {
            log.debug "Missed collect-file cache -- cause: ${e}"
        }
        catch( Exception e ) {
            log.debug "Unable retrieve cached collect-files from: ${temp}", e
            temp.delete()
        }
        return null
    }

    private void cacheCollectedFile(List items, Path temp) {
        def copy = new ArrayList(items)
        try {
            KryoHelper.serialize(copy, temp)
            log.debug "Saved collect-files list to: ${temp}"
        }
        catch( Exception e ) {
            log.warn "Cannot cache collected files -- See the log file for details"
            log.debug("Error serialising collected files", e)
        }

    }

}
