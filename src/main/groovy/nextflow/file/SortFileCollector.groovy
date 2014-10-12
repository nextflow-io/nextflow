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
import net.openhft.chronicle.map.ChronicleMap
import net.openhft.chronicle.map.ChronicleMapBuilder
import nextflow.sort.ChronicleSort
import nextflow.util.KryoHelper

/**
 *  Helper class used to aggregate values having the same key
 *  to files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SortFileCollector extends FileCollector implements Closeable {

    /**
     * Hold a single entry in the tree-set
     */
    static class IndexEntry implements Serializable {

        /*
         * User provided grouping key used to sort entries in the (secondary) index
         */
        Comparable group

        /*
         * The entry index i.e. its key in the store map
         */
        long index

        // required by kryo de-serialization
        protected IndexEntry () {}

        IndexEntry( Comparable a, long n ) {
            this.group = a
            this.index = n
        }

    }

    static class IndexSort implements Comparator<IndexEntry> {

        final SortFileCollector collector;

        final Closure<Comparable> sort

        IndexSort( SortFileCollector obj ) {
            this.collector = obj
            this.sort = obj.sort
        }

        @Override
        int compare(IndexEntry e1, IndexEntry e2) {

            def k = e1.group <=> e2.group
            if( k != 0 )
                return k

            if( sort ) {
                def v1 = collector.getDataAt(e1.index)
                def v2 = collector.getDataAt(e2.index)
                return sort.call(v1) <=> sort.call(v2)
            }

            return e1.index <=> e2.index
        }
    }


    static final OpenOption[] APPEND = [StandardOpenOption.CREATE, StandardOpenOption.APPEND, StandardOpenOption.WRITE] as OpenOption[]

    public Closure sort

    public long entries

    private long count

    private File tempDir;

    private ChronicleMap<Long,Object> store;

    private ChronicleSort<IndexEntry> index;

    /**
     * Creates a file collector
     *
     * @param path The path the temporary files will be stored
     */
    SortFileCollector( Path path = null ) {
        if( !path ) {
            tempDir = Files.createTempDirectory('nxf-collect').toFile()
        }
        else {
            // read the path attributes
            def attr = Files.readAttributes(path, BasicFileAttributes)
            // when it is a dir => create a temp file there
            if ( attr.isDirectory() ) {
                tempDir = path.toFile()
            }
            // if its a file use it
            else if( attr.isRegularFile()) {
                throw new FileAlreadyExistsException("Cannot create sort path: $path -- A path with the same name already exists")
            }
            else {
                tempDir = Files.createDirectories(path).toFile()
            }
        }
    }

    /**
     * Creates the MapDB data structures
     */
    private void createStoreAndIndex() {

        /*
        * Note: on Mac OSX memory mapped files are allocated eagerly, thus using a number of
        * expected entries too big will result in a very big file to be created in the file system
        *
        * See:
        * https://github.com/OpenHFT/Chronicle-Map/blob/master/README.md#size-of-space-reserved-on-disk
        * https://groups.google.com/d/msg/java-chronicle/6UurY1qscS8/GZd28kjSLCsJ
        */
        if( !entries ) {
            entries = "Mac OS X".equals(System.getProperty("os.name")) ? 1_000_000L : 5_000_000_000L;
        }

        File data = new File(tempDir, "store.dat")
        store = ChronicleMapBuilder.of(Long,Object).entries(entries).create(data)

        index = new ChronicleSort<>()
                    .entries(entries)
                    .tempDir(tempDir)
                    .comparator( new IndexSort(this) )
                    .create() as ChronicleSort

    }

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
     * Append a user value to the file collection
     *
     * @param key The grouping key
     * @param value The value to append
     * @return The {@link SortFileCollector} self object
     */
    @Override
    SortFileCollector append( String key, value ) {
        if( !store )
            createStoreAndIndex()

        def bytes = KryoHelper.serialize(value);
        store.put( count, bytes )
        index.add( new IndexEntry(group: key, index: count) )
        count++

        return this
    }


    Object getDataAt( long index ) {
        def bytes = (byte[])store.get(index);
        KryoHelper.deserialize(bytes)
    }


    /**
     * Append the content of a file to the target file having {@code key} as name
     *
     * @param key
     * @param fileToAppend
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

        if( !store ) {
            // nothing  to do
            return
        }

        def last = null
        OutputStream output = null

        index.sort { IndexEntry entry ->
            def name = entry.group
            if( last != name ) {
                // close current output stream
                output?.closeQuietly()
                // set the current as the 'next' last
                last = name

                /*
                 * given a 'group' name returns the target file where to save
                 */
                def target = closure.call(name)
                output = Files.newOutputStream(target, APPEND)

                /*
                 * write the 'seed' value
                 */
                if( seed instanceof Map ) {
                    if( seed.containsKey(name) ) appendStream(normalizeToStream(seed.get(name)), output)
                }
                else if( seed ) {
                    appendStream(normalizeToStream(seed), output)
                }

            }

            /*
             * add the current value
             */
            def bytes = (byte[])store.get(entry.index)
            def val = KryoHelper.deserialize(bytes)
            appendStream(normalizeToStream(val), output)

        }

        // close last output stream
        output?.closeQuietly()
    }

    @Override
    void close() {
        log.debug "Closing collector - begin"
        store?.closeQuietly()
        index?.closeQuietly()
        log.debug "Closing collector - complete"
    }

}
