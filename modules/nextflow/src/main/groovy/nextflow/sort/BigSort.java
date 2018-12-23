/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.sort;

import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import groovy.lang.Closure;
import nextflow.extension.FilesEx;
import nextflow.file.FileHelper;
import nextflow.util.MemoryUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Implementation of External sorting
 * http://en.wikipedia.org/wiki/External_sorting
 * <p>
 * Items to be sorted a stored to off-memory key-value store or splitted in slices
 * having a maximum amount of entries
 * <p>
 * Each slice is ordered in memory and an index file is stored into a temporary file
 * <p>
 * When all slices are sorted they are merged into the final sort
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract public class BigSort<V> implements Closeable {

    private static final Logger log = LoggerFactory.getLogger(BigSort.class);

    /**
     * number of items add to the collection to be sort
     */
    private long count = 0;

    /**
     * Folder where temporary data is stored
     */
    private Path tempDir;

    /**
     * List of index files
     */
    private final List<Path> indexFiles = new LinkedList<>();

    /**
     * List containing the beginning offset of slices after the first
     */
    private final List<Long> slices = new LinkedList<>();

    private long sliceMaxSize = new MemoryUnit("200 MB").toBytes();

    private int sliceMaxItems = 250_000;

    private long sliceSize = 0;

    private boolean deleteTempFilesOnClose = true;

    private Comparator comparator;

    /**
     * @param value Set {@code true} to delete the temporary files created by the storing on termination
     *              or {@code false} to keep on the file system (default: {@code true})
     * @return The {@link BigSort} instance itself
     */
    public BigSort deleteTempFilesOnClose(boolean value) {
        this.deleteTempFilesOnClose = value;
        return this;
    }

    /**
     * @param value Maximum bytes size of a sorting slice
     * @return The {@link BigSort} instance itself
     */
    public BigSort sliceMaxSize( long value ) {
        this.sliceMaxSize = value;
        return this;
    }

    /**
     * @param value Maximum number of items in each sorting slice
     * @return The {@link BigSort} instance itself
     */
    public BigSort sliceMaxItems( int value ) {
        this.sliceMaxItems = value;
        return this;
    }

    /**
     * @param folder Folder where temporary files will be stored. If the directory
     *               does not exist, it is created automatically.
     * @return The {@link BigSort} instance itself
     */
    public BigSort tempDir( Path folder ) {
        tempDir = folder;
        return this;
    }

    /**
     * @return Folder where temporary files will be stored
     */
    public Path getTempDir() {
        return tempDir;
    }

    /**
     * @param c {@link java.util.Comparator} used to sort collection entries
     * @return The {@link BigSort} instance itself
     */
    public BigSort<V> comparator( Comparator c ) {
        this.comparator = c;
        return this;
    }

    public Comparator getComparator() { return comparator; }

    public long getSliceMaxSize() { return sliceMaxSize; }

    public long getSliceMaxItems() { return sliceMaxItems; }

    public boolean getDeleteTempFilesOnClose() { return deleteTempFilesOnClose; }

    List<Long> getSlices() { return Collections.unmodifiableList(slices); }

    /**
     * Creates and initialize temporary data structures required by the sorting
     *
     * @return The {@link BigSort} instance itself
     * @throws IOException
     */
    public BigSort create() throws IOException {
        if( tempDir == null ) {
            tempDir = FileHelper.createLocalDir("bigsort");
        }
        else if ( !Files.exists(tempDir) )  {
            Files.createDirectories(tempDir);
        }

        if( log.isTraceEnabled())
            log.trace("BigSort sliceMaxSize: {}; sliceMaxItems: {}; temp dir: {}",sliceMaxSize, sliceMaxItems, tempDir);

        if( comparator == null ) {
            comparator = DefaultComparator.INSTANCE;
        }

        return this;
    }

    /**
     * Put an entry into the underlying storage
     *
     * @param key   The entry key
     * @param value The entry itself
     * @return The number of bytes required to store this entry
     */
    abstract protected int put(long key, V value);

    /**
     * Retrieve an entry from the underlying storage
     *
     * @param key The entry key
     * @return The entry object
     */
    abstract protected V get(long key);

    /**
     * Add an entry to the collection to be sorted
     *
     * @param value
     */
    public void add(V value) {
        sliceSize += put(count, value);
        count++;

        if (sliceSize > sliceMaxSize || (count % sliceMaxItems) == 0) {
            slices.add(count);
            sliceSize = 0;
        }

        if (log.isTraceEnabled() && (count % 10_000) == 0) {
            log.trace(stats());
        }
    }

    /**
     * Sort the collection of added entries, invoking the user provided {@link groovy.lang.Closure} for each of them
     *
     * @param closure
     */
    public void sort(Closure closure) {

        int index = 0;
        long start = 0;
        try {
            long begin = System.currentTimeMillis();

            /*
             * Sort each slice
             */
            for (long pos : slices) {
                sliceSort(index++, start, pos);
                start = pos;
            }

            // the last chunk
            if (start < count) {
                sliceSort(index++, start, count);
            }

            /*
             * Sort all the partial sorted slices
             */
            long end2 = System.currentTimeMillis();
            double delta1 = (double)(end2 - begin) /1000;        // time required for slices sorting

            externalSort(closure);

            long end3 = System.currentTimeMillis();
            double delta2 = (double)(end3 - end2) /1000;         // time required for external sort
            double delta3 = (double)(end3 - begin) /1000;       // total time

            log.debug("Sort completed -- entries: {}; slices: {}; internal sort time: {} s; external sort time: {} s; total time: {} s", count, slices.size()+1, delta1, delta2, delta3);

        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Given a big number of entries to sort, this method sort a slice of them
     * between the {@code start} and {@code end} indexes
     * <p>
     * The sorted array of indexes are save to a temporary file
     *
     * @param sliceIndex The n-th index of the slice to sort
     * @param start Index from where start the sorting (including it)
     * @param end Index of the last entry to sort (excluding it)
     * @throws IOException
     */
    protected void sliceSort(int sliceIndex, long start, long end) throws IOException {
        assert start < end;

        int len = (int) (end - start);
        Object[] buffer = new Object[len];

        // initialize the array buffer
        for (int i = 0; i < len; i++)
            buffer[i] = i;

        // sort it by using the custom comparator
        Arrays.sort(buffer, new SliceComparator(this, start));

        // add to the list of all index files
        Path file = tempDir.resolve("slice-" + String.valueOf(sliceIndex));
        indexFiles.add(file);

        // save the sorted buffer to a file
        DataOutputStream out = new DataOutputStream(Files.newOutputStream(file));

        // 1. write the offset
        out.writeLong(start);

        // 2. save all the buffer sorted content
        for (int i = 0; i < len; i++)
            out.writeInt((int) buffer[i]);

        // 3. close the file
        out.close();
    }


    /**
     * External sort uses the sorted index created by previous step and
     * merge the final sorting, picking the lowest value from each sorted slice (lane)
     *
     * @param closure
     * @throws IOException
     */
    @SuppressWarnings("unchecked")
    protected void externalSort(Closure closure) throws IOException {
        final int len = indexFiles.size();
        List<Lane> lanes = new ArrayList<>(len);

        /*
         * 1. initialize all lanes
         */
        for (Path file : indexFiles) {
            lanes.add( new Lane(this, new DataInputStream(Files.newInputStream(file))) );
        }

        /*
         * 2. sort by the first value in each lane
         */
        Collections.sort(lanes, new Comparator<Lane>() {
            @Override
            public int compare(Lane o1, Lane o2) {
                return comparator.compare(o1.value, o2.value);
            }

        });

        while (!lanes.isEmpty()) {

            /*
             * 3. pick the lowest value
             */
            Lane cursor = lanes.get(0);

            /*
             * 4. the user provided closure to that it can process the entry found
             */
            closure.call(cursor.value);

            try {
                /*
                 * 5. read next value
                 */
                Object newVal = cursor.next();

                /*
                 * 6. keep the lanes list sorted, moving up the 'newVal' just picked
                 */
                int i=1;
                while( i<lanes.size() && comparator.compare(newVal, lanes.get(i).value) > 0 ) {
                    i++;
                }

                if( i>1 ) {
                    lanes.add(i, cursor);
                    lanes.remove(0);
                }
            }
            catch (EOFException e) {
                /*
                 * 7. remove the 'lane' item when it gets to the end of the index file
                 */
                cursor.close();
                lanes.remove(cursor);
            }
        }

    }

    /**
     * A simple stats information string
     */
    public String stats() {
        return "Added " + String.valueOf(count) + " lines";
    }

    /**
     * Remove all temporary files
     */
    public void close() {
        if( deleteTempFilesOnClose ) {
            FilesEx.deleteDir(tempDir);
        }
    }


    /**
     * Custom comparator to sort the entries in each slice by using the
     * {@link BigSort#comparator} object
     */
    static class SliceComparator implements Comparator {

        private long offset;

        private long count;

        private long cacheHits;

        private long cacheMissed;

        private long begin;

        private final BigSort sort;

        private final Cache<Integer, Object> cache;

        public SliceComparator(BigSort sort, long offset) {
            this.sort = sort;
            this.offset = offset;
            this.cache =  new Cache<>(sort.sliceMaxItems);
            this.begin = System.currentTimeMillis();
        }

        protected Object get(Integer key) {
            Object result = cache.get(key);

            if (result != null) {
                cacheHits++;
            }
            else {
                result = sort.get(offset + key);
                cache.put(key, result);
                cacheMissed++;
            }

            return result;
        }

        /**
         * Implements the compare strategy
         *
         * @param key1
         * @param key2
         * @return
         */
        @Override
        @SuppressWarnings("unchecked")
        public int compare(Object key1, Object key2) {

            int i1 = (Integer) key1;
            int i2 = (Integer) key2;

            Object v1 = get(i1);
            Object v2 = get(i2);

            if( log.isTraceEnabled() ) {
                count++;
                if( (count % 1_000_000) == 0 )
                    log.trace(stats());
            }

            int ret = sort.comparator.compare(v1, v2);
            if (ret == 0) {
                ret = i1 - i2;
            }

            return ret;
        }

        /**
         * @return Some runtime information about the on-going comparison
         */
        private String stats() {
            final long delta = System.currentTimeMillis() - begin;
            begin = System.currentTimeMillis();
            return "Offset: " + String.valueOf(offset) + " - Comparisons " + String.valueOf(count) + " done -- cache hits: " + String.valueOf(cacheHits) + " - missed: " + String.valueOf(cacheMissed) + " -- Elapsed: " + String.valueOf((double) delta / 1000) + " secs";
        }

    }

    /**
     * Implements a simple MRU cache
     *
     * @param <K>
     * @param <V>
     */
    static class Cache<K, V> extends LinkedHashMap<K, V> {

        private int max;

        public Cache(int size) {
            this.max = size;
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > max;
        }

    }

    /**
     * Model a cursor over a slice index file. Every time the method {@link #next()} is invoked
     * a new position entry is read from the index file, the position information plus the offset
     * if used to lookup a value in the main store map
     */
    static class Lane implements Closeable {

        private final DataInputStream stream;
        private final BigSort sort;

        private Object value;
        private long offset;
        private int pos;

        Lane(BigSort sort, DataInputStream stream1) throws IOException {
            this.sort = sort;
            this.stream = stream1;
            this.offset = stream.readLong();
            next();
        }

        public Object next() throws IOException  {
            pos = stream.readInt();
            value = sort.get(offset + pos);
            return value;
        }

        public Object getValue() { return value; }

        @Override
        public void close() throws IOException {
            stream.close();
        }
    }

    /**
     * Default comparator simply assumes that comparing object implements
     * {@link java.lang.Comparable} interface
     */
    static class DefaultComparator implements Comparator {

        static DefaultComparator INSTANCE = new DefaultComparator();

        @Override
        @SuppressWarnings("unchecked")
        public int compare(Object o1, Object o2) {
            return ((Comparable)o1).compareTo(o2);
        }
    }


}