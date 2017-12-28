/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.util

import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.ZipEntry
import java.util.zip.ZipInputStream
import java.util.zip.ZipOutputStream

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.FileHelper

/**
 * Make a script executable in a remote session, embedding the script classes and external libraries.
 * <p>
 * Given a session it traverse and compress all the files in the paths by the session property
 * {@link Session#libDir}.
 * <p>
 * The compressed binary is serialized and transported to the remote node where they can be un-compressed
 * and made available on the remote classpath
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class RemoteSession implements Serializable, Closeable {

    private static final long serialVersionUID = 8675142190223473251L ;

    final UUID sessionId

    final List<byte[]> libs = []

    final boolean statsEnabled

    transient boolean isLibInitialized

    @Lazy
    transient List<Path> libDir = lazyLibDir()

    @Lazy
    transient List<Path> classpath = { ConfigHelper.resolveClassPaths(libDir) }()

    protected RemoteSession() { }

    /**
     * Creates the object and compress the classes and libraries ina binary form
     *
     * @param session A {@link Session} object for which the {@link Session#libDir} has been initialized
     */
    RemoteSession(Session session) {
        sessionId = session.getUniqueId()
        statsEnabled = session.statsEnabled
        for( Path entry : session.getLibDir() ) {
            libs.add( zip(entry) )
        }
    }

    /**
     * @param path A {@link File} representing either a regular file or a directory
     * @return The source path as a string. It guarantees that a directory path ends with a slash character
     */
    private String pathEndWithSeparator( Path path ) {
        def result = path.complete().toString()

        if( !path.isDirectory() )
            return result

        if( !result.endsWith(File.separator))
            result += File.separator

        result
    }

    /**
     * Zip the folder content recursively and return a byte[] object
     *
     * @param dir The {@link File} folder to zip
     * @return A byte array holding the folder content as a compressed binary array
     */
    protected byte[] zip( Path dir ) {
        assert dir
        assert dir.isDirectory()

        int count = 0
        def base = pathEndWithSeparator(dir)
        def buffer = new ByteArrayOutputStream()
        ZipOutputStream zip = new ZipOutputStream(buffer);

        dir.eachFileRecurse { Path file ->

            String name = pathEndWithSeparator(file)
            if( name.startsWith(base)) {
                name = name.substring(base.length());
            }
            ZipEntry entry = new ZipEntry(name);
            zip.putNextEntry(entry);

            if( file.isFile() ) {
                count++
                // append the file content
                Files.copy(file, zip)
            }
            // Complete the entry
            zip.closeEntry();
        }

        zip.close();
        buffer.close()
        log.debug "Packed library ($count files) at path: $dir"
        return buffer.toByteArray()
    }

    /**
     * Unzip a byte array buffer into a target {@link File} folder
     *
     * @param bytes A byte array
     * @param target
     * @return
     */
    protected Path unzip( byte[] bytes, Path target = null ) {
        assert bytes

        if( !target )
            target = FileHelper.createLocalDir()

        int count=0
        ZipInputStream zip = new ZipInputStream(new ByteArrayInputStream(bytes))
        ZipEntry entry
        while( (entry=zip.getNextEntry()) != null ) {

            def file = target.resolve(entry.getName());
            if(entry.isDirectory()) {
                continue
            }

            if( !file.parent.exists() && !file.parent.mkdirs() )
                throw new IllegalStateException("Unable to create target directory: ${file.parent} -- check file permissions")

            count++
            Files.copy(zip, file)
        }
        zip.close();
        log.debug "Staged library ($count files) to path: $target"

        return target
    }

    /**
     * @return The list of folder containing the unzip libraries
     */
    protected List<Path> lazyLibDir() {
        def result = (List<Path>)libs.collect { byte[] it -> unzip(it) }
        isLibInitialized = true
        return result
    }

    /**
     * Close the session i.e. delete all temporary directories
     *
     * @throws IOException
     */
    @Override
    void close() throws IOException {
        if(!isLibInitialized) return
        libDir.each { Path it -> it.deleteDir() }
    }
}
