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
import java.nio.file.Files
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.StandardOpenOption
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
/**
 *  Helper class used to aggregate values having the same key
 *  to files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileCollector implements Closeable {

    static final OpenOption[] APPEND = [StandardOpenOption.APPEND, StandardOpenOption.WRITE] as OpenOption[]

    public Boolean newLine

    public seed

    protected Path temp

    protected ConcurrentMap<String,Path> cache = new ConcurrentHashMap<>()

    FileCollector( Path store = null ) {
        if( !store )
            store = Files.createTempDirectory('nxf-clt')

        else {
            store.createDirIfNotExists()
        }
        this.temp = store
    }

    protected Path _file( String name ) {
        (Path)cache.getOrCreate(name) {
            def result = Files.createFile(temp.resolve(name))
            if( seed instanceof Map && seed.containsKey(name)) {
                append0(result, _value(seed.get(name)))
            }
            else if( seed ) {
                append0(result, _value(seed))
            }
            return result
        }
    }

    protected InputStream _value( value ) {
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


    FileCollector append( String key, value ) {
        append0( _file(key), _value(value))
        return this
    }


    /**
     * Append the content of a file to the target file having {@code key} as name
     *
     * @param key
     * @param fileToAppend
     */
    protected void append0( Path source, InputStream stream ) {
        int n
        byte[] buffer = new byte[10 * 1024]
        def output = Files.newOutputStream(source, APPEND)

        try {
            while( (n=stream.read(buffer)) > 0 ) {
                output.write(buffer,0,n)
            }
            // append the new line separator
            if( newLine )
                output.write( System.lineSeparator().bytes )
        }
        finally {
            stream.closeQuietly()
            output.closeQuietly()
        }
    }

    /**
     *
     * @return The number of files in the appender accumulator
     */
    int size() {
        cache.size()
    }

    boolean isEmpty() {
        cache.isEmpty()
    }

    boolean containsKey(String key) {
        return cache.containsKey(key)
    }

    Path get(String name) {
        cache.get(name)
    }

    List<Path> getFiles() {
        new ArrayList<Path>(cache.values())
    }

    List<Path> moveFiles(Path target) {
        target.createDirIfNotExists()

        def result = []
        moveFiles { String name ->
            Path newFile = target.resolve(name)
            result << newFile
            return newFile
        }
        return result
    }

    void moveFiles( Closure<Path> closure ) {

        def result = []
        Iterator<Path> itr = cache.values().iterator()
        while( itr.hasNext() ) {
            def item = itr.next()
            def target = closure.call(item.getName())
            result << Files.move(item, target)
            itr.remove()
        }

    }

    @Override
    void close() {
        temp?.deleteDir()
    }


}
