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
import static org.iq80.leveldb.impl.Iq80DBFactory.factory

import java.nio.ByteBuffer
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.sort.LevelDbSort
import nextflow.util.KryoHelper
import org.iq80.leveldb.DB
import org.iq80.leveldb.Options
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

        final DB store;

        final Closure<Comparable> sort

        IndexSort( DB store, Closure<Comparable> sort ) {
            this.store = store
            this.sort = sort
        }

        @Override
        int compare(IndexEntry e1, IndexEntry e2) {

            def k = e1.group <=> e2.group
            if( k != 0 )
                return k

            if( this.sort ) {
                def v1 = getValue(e1.index)
                def v2 = getValue(e2.index)
                return sort.call(v1) <=> sort.call(v2)
            }

            return e1.index <=> e2.index
        }

        Object getValue( long index ) {
            def raw = (byte[])store.get(bytes(index));
            KryoHelper.deserialize(raw)
        }

    }


    Closure sort

    Long sliceMaxSize

    Integer sliceMaxItems

    private long count

    private DB store;

    private LevelDbSort<IndexEntry> index;

    private Path fTempDir

    /**
     * Creates the underlying map data structures
     */
    private void createStoreAndIndex() {
        log.trace "Creating sort file collector -- temp dir: $tempDir"

        def storeDir = getTempDir().resolve("data").toFile()
        Options options = new Options().createIfMissing(true);
        store = factory.open(storeDir, options);

        def indexDir = getTempDir().resolve('index')
        def result = new LevelDbSort<IndexEntry>()
        result.comparator( new IndexSort(store, this.sort) )
        result.tempDir(indexDir)
        result.deleteTempFilesOnClose(this.deleteTempFilesOnClose)

        if( sliceMaxSize ) result.sliceMaxSize(sliceMaxSize)
        if( sliceMaxItems ) result.sliceMaxItems(sliceMaxItems)

        index = result.create()

        fTempDir = getTempDir()
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
    SortFileCollector add( String key, value ) {

        // allocate data structures
        if( store == null ) { createStoreAndIndex() }

        // serialise the main value
        def payload = KryoHelper.serialize(value);

        store.put( bytes(count), payload )
        index.add( new IndexEntry(group: key, index: count) )

        count++

        return this
    }

    private static byte[] bytes(long value) {
        return ByteBuffer.allocate(8).putLong(value).array();
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

    /**
     * {@inheritDoc}
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
     * {@inheritDoc}
     */
    void saveFile( Closure<Path> closure ) {

        if( !index )
            return

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
            def bytes = (byte[])store.get(bytes(entry.index))
            def val = KryoHelper.deserialize(bytes)
            appendStream(normalizeToStream(val), output)
        }

        // close last output stream
        output?.closeQuietly()
    }

    /**
     * Close sorting structures instance and cleanup temporary files
     */
    @Override
    void close() {
        log.trace "Closing sorting dbs"
        store?.closeQuietly()
        index?.closeQuietly()

        // finally invoke the the parent close
        super.close()

        if( deleteTempFilesOnClose && fTempDir ) {
            fTempDir.deleteOnExit()
        }
        else {
            log.debug "FileCollector temp dir not removed: $tempDir"
        }
    }

}
