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
import static org.iq80.leveldb.impl.Iq80DBFactory.factory

import java.nio.ByteBuffer
import java.nio.file.Files
import java.nio.file.Path

import com.google.common.hash.HashCode
import com.google.common.hash.Hasher
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.sort.LevelDbSort
import nextflow.util.CacheHelper
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

        /**
         * The entry hash code
         */
        byte[] hash

        // required by kryo de-serialization
        protected IndexEntry () {}

        IndexEntry( Comparable a, long n, byte[] h ) {
            this.group = a
            this.index = n
            this.hash = h
        }

    }

    static class IndexSort implements Comparator<IndexEntry> {

        final DB store;

        final Closure<Comparable> sort

        final Comparator comparator

        IndexSort( DB store ) {
            this.store = store
            this.sort = null
            this.comparator = null
        }

        IndexSort( DB store, Closure<Comparable> criteria ) {
            this.store = store
            this.sort = criteria
            this.comparator = null
        }

        IndexSort( DB store, Comparator criteria ) {
            this.store = store
            this.comparator = criteria
            this.sort = null
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

            if( this.comparator ) {
                def v1 = getValue(e1.index)
                def v2 = getValue(e2.index)
                return comparator.compare(v1,v2)
            }

            return e1.index <=> e2.index
        }

        Object getValue( long index ) {
            def raw = (byte[])store.get(bytes(index));
            KryoHelper.deserialize(raw)
        }

    }

    def sort

    Long sliceMaxSize

    Integer sliceMaxItems

    private long count

    private DB store;

    private LevelDbSort<IndexEntry> index

    private Path tempDir

    private SequentialFileStore trail


    void setSort( def value ) {
        if( value == null || value instanceof Closure || value instanceof Comparator )
            this.sort = value
        else
            throw new AbortOperationException("Not a valid `sort` criteria: $sort -- It can be a Closure or a Comparator object")
    }

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
        result.comparator(createSortComparator())
        result.tempDir(indexDir)
        result.deleteTempFilesOnClose(this.deleteTempFilesOnClose)

        if( sliceMaxSize ) result.sliceMaxSize(sliceMaxSize)
        if( sliceMaxItems ) result.sliceMaxItems(sliceMaxItems)

        index = result.create()
        tempDir = getTempDir()

    }

    private IndexSort createSortComparator() {
        if( !sort )
            return new IndexSort(store)

        if( sort instanceof Closure )
            return new IndexSort(store, this.sort as Closure)

        if( sort instanceof Comparator )
            return new IndexSort(store, this.sort as Comparator)

        throw new IllegalStateException("Not a valid `sort` criteria object: $sort")
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
        def hash = hash(key, value)

        store.put( bytes(count), payload )
        index.add( new IndexEntry(key, count, hash) )

        count++

        return this
    }

    private static byte[] bytes(long value) {
        return ByteBuffer.allocate(8).putLong(value).array();
    }

    private byte[] hash(String key, value) {
        def h
        h = CacheHelper.hasher(key, hashMode)
        h = CacheHelper.hasher(h, value, hashMode)
        return h.hash().asBytes()
    }

    @Override
    protected HashCode makeHash() {

        if( !index )
            return null

        Hasher hasher = cacheable ? CacheHelper.hasher( hashKeys, CacheHelper.HashMode.STANDARD ) : null
        if( hasher ) {
            log.trace "  hasher: ${CacheHelper.hasher( hashKeys, CacheHelper.HashMode.STANDARD ) } \n"
        }

        // index route file
        trail = new SequentialFileStore(tempDir.resolve('sort.index'))

        /*
         * Sort the index and record the entry indices into the 'trail' structure.
         */
        def last = null
        index.sort { IndexEntry entry ->
            // the index is sorted by group, when the current 'group'
            // is different from the previous one save it to the 'trail'
            if( last != entry.group ) {
                last = entry.group
                def bytes = KryoHelper.serialize(last)
                trail.writeLong( -bytes.length )
                trail.writeBytes( bytes )
            }
            // save the entry position on the index
            trail.writeLong(entry.index)

            if( hasher ) {
                hasher = CacheHelper.hasher(hasher, entry.hash, CacheHelper.HashMode.STANDARD)
                log.trace "  index: $entry.index - ${CacheHelper.hasher(entry.hash).hash()} \n"
            }

        }

        return hasher?.hash()
    }

    /**
     * {@inheritDoc}
     */
    void saveFile( Closure<Path> closure ) {

        if( !trail )
            return

        /*
         * turn the index buffer into 'read' mode
         */
        trail.flip()

        long index = -1
        def last = null
        def name = null
        OutputStream output = null

        /*
         * the 'trail' contain the ordered sequence of item indices
         */
        while( true ) {
            /*
             * read the next item index
             */
            try {
                index = trail.readLong()
            }
            catch( EOFException e ) {
                break
            }

            /*
             * when it's a negative value the next element is the grouping element
             * and the value represent the number of bytes to read
             */
            if( index < 0  ) {
                def data = trail.readBytes((int)-index)
                name = KryoHelper.deserialize(data)
                continue
            }

            /*
             * fetch the next value from the store (this also may retried the header)
             */
            def bytes = (byte[])store.get(bytes(index))
            def val = KryoHelper.deserialize(bytes)
            def data = normalizeToStream(val)

            /*
             * check if a new group has been found
             */
            if( last != name ) {
                // close current output stream
                output?.closeQuietly()
                // set the current as the 'next' last
                last = name

                /*
                 * given a 'group' name returns the target file where to save
                 */
                def target = closure.call(name)
                output = Files.newOutputStream(target, TRUNCATE)

                /*
                 * append the the header
                 */
                appendHeader(data, name, output)
            }

            /*
             * append the current item value
             */
            appendStream(data, output)
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
        trail?.closeQuietly()

        // finally invoke the the parent close
        super.close()

        if( deleteTempFilesOnClose && tempDir ) {
            tempDir.deleteDir()
        }
        else {
            log.debug "FileCollector temp dir not removed: $tempDir"
        }
    }

}
