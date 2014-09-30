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

package test
import static org.iq80.leveldb.impl.Iq80DBFactory.factory

import java.nio.ByteBuffer

import groovy.transform.CompileStatic
import nextflow.util.KryoHelper
import nextflow.util.MemoryUnit
import org.iq80.leveldb.DB
import org.iq80.leveldb.Options

@CompileStatic
class TestLevelDb5 {

    static byte[] bytes(int value) {
        ByteBuffer.allocate(4).putInt(value).array();
    }

    int count = 0
    final DB data
    final File indexFolder
    final List<File> indexFiles = new LinkedList<>()
    List<Long> slices = new LinkedList<>()

    long maxSize = new MemoryUnit('100 MB').toBytes()
    int maxItems = 150_000
    long size = 0

   String asString(byte[] bytes) {
        def result = new StringBuilder()
        for (byte b : bytes) {
            result << String.format("0x%x ", b)
        }
        result.toString()
    }


    def TestLevelDb5() {

        Options options = new Options()
                .createIfMissing(true)
        data = factory.open(new File("/Users/pditommaso/Downloads/db/data"), options);

        indexFolder = new File("/Users/pditommaso/Downloads/db/index")
        indexFolder.mkdirs()
    }

    def add( String line ) {
        //println "Adding $count"
        def key = bytes(count)
        def value = KryoHelper.serialize(line)
        data.put(key, value )
        count++

        // increase size count
        size += 4
        size += value.size()

        if( size > maxSize || (count % maxItems) == 0 ) {
            slices.add( (long)count )
            size = 0
        }

        if( count % 10_000 == 0 )
            println stats()
    }

    def void sort(Closure closure) {
        int index = 0
        long start = 0

        for( long pos : slices ) {
            sliceSort(index++, start, pos)
            start = pos
        }

        // the last chunk
        if( start < count ) {
            sliceSort(index++, start, count)
        }

        externalSort(closure)
    }

    def void sliceSort( int index, long start, long end ) {

        assert start < end
        def len = end-start
        Object[] buffer = new Object[len]
        for( int i=0; i<len; i++ )
            buffer[i] = new Integer(i)

        Arrays.sort(buffer, new MyComparator((int)start))

        def file = new File(indexFolder, "ndx-$index")
        indexFiles << file

        def out = new DataOutputStream(new FileOutputStream(file))
        // write the offset
        out.writeInt((int)start)

        for( int i=0; i<len; i++ )
            out.writeInt( buffer[i] as Integer )

        out.close()

    }

    def void externalSort( Closure closure ) {
        long begin = System.currentTimeMillis()
        println "Starting external sort -- ${indexFiles.size()} entries"
        DataInputStream[] streams = new DataInputStream[indexFiles.size()]
        int[] indexes = new int[indexFiles.size()]
        int[] offsets = new int[indexFiles.size()]
        Comparable[] values = new Comparable[indexFiles.size()]

        int i=0
        for( File file : indexFiles ) {
            streams[i] = new DataInputStream( new FileInputStream(file) )
            offsets[i] = streams[i].readInt()
            indexes[i] = streams[i].readInt()
            values[i] = (Comparable)KryoHelper.deserialize(data.get( bytes(offsets[i] + indexes[i]) ))
            i++
        }


        int p = -1
        def len = streams.length
        def terminated = 0
        while( terminated<len ) {
            Comparable min = null
            for( i=0; i<values.length; i++ ) {
                if( values[i] == null )
                    continue

                if( min == null ) {
                    min = values[i]
                    p = i
                }

                else if( (min <=> values[i]) > 0 ) {
                    min = values[i]
                    p = i
                }
            }

            if( min == null )
                break

            // call
            closure.call(min)

            // clear current
            min = null

            // read next value
            try {
                def index = offsets[p] + streams[p].readInt()
                values[p] = (Comparable)KryoHelper.deserialize(data.get(bytes(index)))
            }
            catch( EOFException e ) {
                values[p] = null
                streams[p].close()
                streams[p] = null
                terminated ++
            }
        }

        println "Ending external sort -- time ${(System.currentTimeMillis()-begin)/1000} secs"
    }


    class MyComparator implements Comparator {

        int offset
        long num
        long cacheHits
        long cacheMissed
        long begin

        MyComparator(int offset) {
            this.offset = offset
            begin = System.currentTimeMillis()
        }

        final Cache<Integer,Comparable> cache = new Cache<>(maxItems)

        private Comparable get( Integer key ) {
            def result = cache.get(key)
            if( result != null ) {
                cacheHits++
            }
            else {
                def dbKey = bytes(offset+key)
                result = (String)KryoHelper.deserialize(data.get(dbKey))
                cache.put(key,result)
                cacheMissed++
            }

            return result
        }


        @Override
        int compare(def key1, def key2) {

            int i1 = (Integer)key1
            int i2 = (Integer)key2

            def v1 = get(i1)
            def v2 = get(i2)

            //
            num++
            if( num % 1_000_000 == 0 )
                println stats()

            def ret = v1 <=> v2
            if( ret == 0 ) {
                ret = i1 <=> i2
            }
            return ret
        }

        def stats() {
            def delta = System.currentTimeMillis() - begin
            begin = System.currentTimeMillis()
            "Offset: $offset - Comparisons $num done -- cache hits: $cacheHits - missed: $cacheMissed -- Elapsed: ${delta/1000} secs"
        }

    }

    @CompileStatic
    class Cache<K,V> extends LinkedHashMap<K,V> {

        int max

        Cache(int size) {
            this.max = size
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > max
        }

    }

    def stats() {
        "Added $count lines"
    }

    def close() {
        data.close()
    }


    static void main(String[] args) {

        def dbFolder = new File('/Users/pditommaso/Downloads/db/')
        dbFolder.deleteDir()

        def collector = new TestLevelDb5()

        println "Begin"
        def now = System.currentTimeMillis()
        def text = new File('/Users/pditommaso/Downloads/hs_ref_GRCh38_chr1.fa') // hs_ref_GRCh38_chr1.fa
        text.eachLine { String it ->
            collector.add(it)
        }
        println "End adding - time: ${(System.currentTimeMillis()-now)/1000} secs -- Stats: ${collector.stats()}"

        def sortFile = new PrintWriter(new FileWriter(new File(dbFolder,'sort.txt') ))
        collector.sort {
            sortFile.println(it)
        }
        sortFile.close()
        println "End sort - time: ${(System.currentTimeMillis()-now)/1000} secs"

        collector.close();
    }

}