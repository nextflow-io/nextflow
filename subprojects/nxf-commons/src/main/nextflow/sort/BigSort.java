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

package nextflow.sort;

import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import groovy.lang.Closure;
import nextflow.extension.FileEx;
import nextflow.util.MemoryUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract public class BigSort<V> implements Closeable {

    private static final Logger log = LoggerFactory.getLogger(BigSort.class);

    private int count = 0;

    private File tempDir;

    private File indexFolder;

    private final List<File> indexFiles = new LinkedList<>();

    private final List<Long> slices = new LinkedList<>();

    private long sliceMaxSize = new MemoryUnit("100 MB").toBytes();

    private int sliceMaxItems = 150_000;

    private long size = 0;

    private boolean deleteTempFilesOnClose = true;

    public BigSort deleteTempFilesOnClose(boolean value) {
        this.deleteTempFilesOnClose = true;
        return this;
    }

    public BigSort sliceMaxSize( long value ) {
        this.sliceMaxSize = value;
        return this;
    }

    public BigSort sliceMaxItems( int value ) {
        this.sliceMaxItems = value;
        return this;
    }

    public BigSort tempDir( File folder ) {
        tempDir = folder;
        return this;
    }

    public File getTempDir() {
        return tempDir;
    }

    public BigSort create() throws IOException {
        if( tempDir == null ) {
            tempDir = Files.createTempDirectory("bigsort").toFile();
        }
        else if ( !tempDir.mkdirs() )  {
            throw new IOException("Unable to create sort temporary folder: " + tempDir);
        }

        indexFolder = new File(tempDir, "index");
        indexFolder.mkdir();

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
        size += put(count, value);
        count++;

        if (size > sliceMaxSize || (count % sliceMaxItems) == 0) {
            slices.add((long) count);
            size = 0;
        }


        if (log.isDebugEnabled() && (count % 10_000) == 0) {
            log.debug(stats());
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
            externalSort(closure);
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Sort a "slice" of entries
     * @param index The n-th index of the slice be merged
     * @param start The
     * @param end
     * @throws IOException
     */
    protected void sliceSort(int index, long start, long end) throws IOException {
        assert start < end;

        int len = (int) (end - start);
        Object[] buffer = new Object[len];

        // initialize the array buffer
        for (int i = 0; i < len; i++)
            buffer[i] = i;

        // sort it by using the custom comparator
        Arrays.sort(buffer, new SortComparator(this, start));

        // add to the list of all index files
        File file = new File(indexFolder, "slice-" + String.valueOf(index));
        indexFiles.add(file);

        // save the sorted buffer to a file
        DataOutputStream out = new DataOutputStream(new FileOutputStream(file));

        // 1. write the offset
        out.writeInt((int) start);

        // 2. save all the buffer sorted content
        for (int i = 0; i < len; i++)
            out.writeInt((int) buffer[i]);

        // 3. close the file
        out.close();
    }

    /**
     *
     * @param closure
     * @throws IOException
     */
    protected void externalSort(Closure closure) throws IOException {
        final long begin = System.currentTimeMillis();
        final int len = indexFiles.size();

        if (log.isDebugEnabled()) {
            log.debug("Starting external sort -- " + String.valueOf(len) + " entries");
        }

        int[] indexes = new int[len];
        int[] offsets = new int[len];
        DataInputStream[] streams = new DataInputStream[len];
        Comparable[] values = new Comparable[len];

        int i = 0;
        for (File file : indexFiles) {
            streams[i] = new DataInputStream(new FileInputStream(file));
            offsets[i] = streams[i].readInt();
            indexes[i] = streams[i].readInt();
            values[i] = (Comparable) get(offsets[i] + indexes[i]);
            i++;
        }


        int p = -1;

        Integer terminated = 0;
        while (terminated < len) {
            Comparable min = null;
            for (i = 0; i < values.length; i++) {
                if (values[i] == null)
                    continue;

                if (min == null) {
                    min = values[i];
                    p = i;
                } else if ((min.compareTo(values[i])) > 0) {
                    min = values[i];
                    p = i;
                }
            }

            // nothing to do
            if (min == null)
                break;

            /*
             * the user provided closure to that it can process the entry found
             */
            closure.call(min);

            // read next value
            try {
                Integer index = offsets[p] + streams[p].readInt();
                values[p] = (Comparable) get(index);
            }
            catch (EOFException e) {
                streams[p].close();
                streams[p] = null;
                values[p] = null;
                terminated++;
            }
        }

        if( log.isDebugEnabled() ) {
            log.debug("Ending external sort -- time " + String.valueOf((double) (System.currentTimeMillis() - begin) / 1000) + " secs");
        }
    }

    public String stats() {
        return "Added " + String.valueOf(count) + " lines";
    }

    public void close() {
        if( deleteTempFilesOnClose ) {
            FileEx.deleteDir(tempDir.toPath());
        }
    }


    /**
     * Custom comparator to sort the entries in each slice
     */
    static class SortComparator implements Comparator {

        private long offset;

        private long count;

        private long cacheHits;

        private long cacheMissed;

        private long begin;

        private final BigSort owner;

        private final Cache<Integer, Comparable> cache;

        public SortComparator(BigSort owner, long offset) {
            this.owner = owner;
            this.offset = offset;
            this.cache =  new Cache<>(owner.sliceMaxItems);
            this.begin = System.currentTimeMillis();
        }

        protected Comparable get(Integer key) {
            Comparable result = cache.get(key);

            if (result != null) {
                cacheHits++;
            }
            else {
                result = (Comparable)owner.get(offset + key);
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
        public int compare(Object key1, Object key2) {

            int i1 = (Integer) key1;
            int i2 = (Integer) key2;

            Comparable v1 = get(i1);
            Comparable v2 = get(i2);

            if( log.isDebugEnabled() ) {
                count++;
                if( (count % 1_000_000) == 0 )
                    log.debug(stats());
            }

            Integer ret = v1.compareTo(v2);
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

}