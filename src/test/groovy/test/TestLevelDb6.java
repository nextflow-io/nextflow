package test;

import static org.iq80.leveldb.impl.Iq80DBFactory.factory;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import groovy.lang.Closure;
import nextflow.extension.FileEx;
import nextflow.util.KryoHelper;
import nextflow.util.MemoryUnit;
import org.iq80.leveldb.DB;
import org.iq80.leveldb.Options;

public class TestLevelDb6 {

    public static byte[] bytes(int value) {
        return ByteBuffer.allocate(4).putInt(value).array();
    }

    public String asString(Byte[] bytes) {
        StringBuilder result = new StringBuilder();
        for (byte b : bytes) {
            result.append(String.format("0x%x ", b));
        }

        return result.toString();
    }

    public TestLevelDb6() throws IOException {

        Options options = new Options().createIfMissing(true);
        data = factory.open(new File("/Users/pditommaso/Downloads/db/data"), options);

        indexFolder = new File("/Users/pditommaso/Downloads/db/index");
        indexFolder.mkdirs();
    }

    public void add(String line) {
        //println "Adding $count"
        byte[] key = bytes(count);
        byte[] value = KryoHelper.serialize(line);
        data.put(key, value);
        count++;

        // increase size count
        size += 4;
        size += value.length;

        if (size > maxSize || (count % maxItems) == 0) {
            slices.add((long) count);
            size = 0;
        }


        if( (count % 10_000) == 0)
            System.out.println(stats());
    }

    public void sort(Closure closure) {
        int index = 0;
        long start = 0;

        try {
            for (long pos : slices) {
                sliceSort(index++, start, pos);
                start = pos;
            }


            // the last chunk
            if (start < count) {
                sliceSort(index++, start, count);
            }


            externalSort(closure);
        }
        catch( IOException e ) {
            throw new RuntimeException(e);
        }

    }

    public void sliceSort(int index, long start, long end) throws IOException {

        assert start < end;
        int len = (int)(end - start);
        Object[] buffer = new Object[len];
        for (int i = 0; i < len; i++)
            buffer[i] = new Integer(i);

        Arrays.sort(buffer, new MyComparator((int) start));

        File file = new File(indexFolder, "ndx-" + String.valueOf(index));
        indexFiles.add(file);

        DataOutputStream out = new DataOutputStream(new FileOutputStream(file));

        // write the offset
        out.writeInt((int) start);

        for (int i = 0; i < len; i++)
            out.writeInt((int)buffer[i]);

        out.close();

    }

    public void externalSort(Closure closure) throws IOException {
        final long begin = System.currentTimeMillis();
        System.out.println("Starting external sort -- " + String.valueOf(indexFiles.size()) + " entries");
        DataInputStream[] streams = new DataInputStream[indexFiles.size()];
        Integer[] indexes = new Integer[indexFiles.size()];
        Integer[] offsets = new Integer[indexFiles.size()];
        Comparable[] values = new Comparable[indexFiles.size()];

        int i = 0;
        for (File file : indexFiles) {
            streams[i] = new DataInputStream(new FileInputStream(file));
            offsets[i] = streams[i].readInt();
            indexes[i] = streams[i].readInt();
            values[i] = (Comparable) KryoHelper.deserialize(data.get(bytes(offsets[i] + indexes[i])));
            i++;
        }


        int p = -1;
        Integer len = streams.length;
        Integer terminated = 0;
        while (terminated < len) {
            Comparable min = null;
            for (i = 0; i < values.length ;i++){
                if (values[i] == null) continue;

                if (min == null) {
                    min = values[i];
                    p = i;
                } else if ((min.compareTo(values[i])) > 0) {
                    min = values[i];
                    p = i;
                }

            }


            if (min == null) break;

            // call
            closure.call(min);

            // clear current
            min = null;

            // read next value
            try {
                Integer index = offsets[p] + streams[p].readInt();
                values[p] = (Comparable) KryoHelper.deserialize(data.get(bytes(index)));
            }
            catch (EOFException e) {
                values[p] = null;
                streams[p].close();
                streams[p] = null;
                terminated++;
            }

        }

        System.out.println("Ending external sort -- time " + String.valueOf((double) (System.currentTimeMillis() - begin) / 1000) + " secs");
    }

    public String stats() {
        return "Added " + String.valueOf(count) + " lines";
    }

    public void close() {
        FileEx.closeQuietly(data);
    }


    private int count = 0;
    private final DB data;
    private final File indexFolder;
    private final List<File> indexFiles = new LinkedList<File>();
    private List<Long> slices = new LinkedList<Long>();
    private long maxSize = new MemoryUnit("100 MB").toBytes();
    private int maxItems = 150_000;
    private long size = 0;

    public class MyComparator implements Comparator {
        public MyComparator(int offset) {
            this.offset = offset;
            begin = System.currentTimeMillis();
        }

        private Comparable get(Integer key) {
            Comparable result = cache.get(key);
            if (result != null) {
                cacheHits++;
            } else {
                byte[] dbKey = bytes(offset + key);
                result = (String) KryoHelper.deserialize(data.get(dbKey));
                cache.put(key, result);
                cacheMissed++;
            }


            return ((Comparable) (result));
        }

        @Override
        public int compare(Object key1, Object key2) {

            int i1 = (Integer) key1;
            int i2 = (Integer) key2;

            Comparable v1 = get(i1);
            Comparable v2 = get(i2);

            //
            num++;
            if( (num % 1_000_000) == 0)
                System.out.println(stats());

            Integer ret = v1.compareTo(v2);
            if (ret == 0) {
                ret = i1 - i2;
            }

            return ((int) (ret));
        }

        public String stats() {
            final Long delta = System.currentTimeMillis() - begin;
            begin = System.currentTimeMillis();
            return "Offset: " + String.valueOf(offset) + " - Comparisons " + String.valueOf(num) + " done -- cache hits: " + String.valueOf(cacheHits) + " - missed: " + String.valueOf(cacheMissed) + " -- Elapsed: " + String.valueOf((double) delta / 1000) + " secs";
        }

        private int offset;

        private long num;

        private long cacheHits;

        private long cacheMissed;

        private long begin;

        private final Cache<Integer, Comparable> cache = new Cache<Integer, Comparable>(maxItems);
    }

    public class Cache<K, V> extends LinkedHashMap<K, V> {

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
