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

package nextflow.file
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx

/**
 * File collector base class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
abstract class FileCollector implements Closeable {

    static final public OpenOption[] APPEND = [StandardOpenOption.CREATE, StandardOpenOption.APPEND, StandardOpenOption.WRITE] as OpenOption[]

    public Boolean newLine

    public seed

    public boolean deleteTempFilesOnClose = true;

    private Path tempDir

    private boolean created

    private boolean closed;

    protected FileCollector() { }

    /**
     * Create the temporary directory
     * @return A {@code Path} instance to the directory location
     */
    protected Path createTempDir() {

        if( !tempDir ) {
            return tempDir = Files.createTempDirectory('nxf-collect')
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
     * @param closure Closure mapping a group name to the respective groupping file e.g. { name -> Paths.get(name) }
     */
    abstract void saveFile( Closure<Path> closure );

    /**
     * Normalize values to a {@link InputStream}
     *
     * @param value The user provided value
     * @return An {@link InputStream} referring the value
     */
    protected InputStream normalizeToStream( value ) {
        if( value instanceof Path )
            return value.newInputStream()

        if( value instanceof File )
            return value.newInputStream()

        if( value instanceof CharSequence )
            return new ByteArrayInputStream(value.toString().getBytes())

        if( value instanceof byte[] )
            return new ByteArrayInputStream((byte[])value)

        throw new IllegalArgumentException("Not a valid file collector argument [${value.class.name}]: $value")
    }

    /**
     * Save the entries collected grouping them into files whose name is given by the
     * correspondent group key
     *
     * @param target A {@link Path} to the folder where files will be saved. If the folder does not exists, is created
     *              automatically
     * @return The list of files where entries have been saved.
     */
    List<Path> saveTo(Path target) {
        target.createDirIfNotExists()

        def result = []
        saveFile { String name ->
            Path newFile = target.resolve(name)
            result << newFile
            return newFile
        }
        return result
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


}
