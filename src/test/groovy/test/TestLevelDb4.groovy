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

import gnu.trove.map.hash.TCustomHashMap
import gnu.trove.strategy.HashingStrategy
import groovy.transform.CompileStatic
import nextflow.util.KryoHelper
import org.iq80.leveldb.DB
import org.iq80.leveldb.DBComparator
import org.iq80.leveldb.DBIterator
import org.iq80.leveldb.Options

@CompileStatic
class TestLevelDb4 {

    static byte[] bytes(int value) {
        ByteBuffer.allocate(4).putInt(value).array();
    }

    int count = 0
    final byte[] EMPTY = [] as byte[]
    final DB data
    final DB index
    final Cache<String> cache
    long cacheHits
    long cacheMissed

   String asString(byte[] bytes) {
        def result = new StringBuilder()
        for (byte b : bytes) {
            result << String.format("0x%x ", b)
        }
        result.toString()
    }


    class MyComparator implements DBComparator{

        @Override
        String name() {
            return 'MyComp'
        }

        private Comparable get( byte[] key ) {
            def result = (Comparable)cache.get(key)
            if( result != null ) {
                cacheHits++
            }
            else {
                result = (String)KryoHelper.deserialize(data.get(key))
                cache.put(key,result)
                cacheMissed++
            }

            return result
        }


        @Override
        int compare(byte[] key1, byte[] key2) {
            def v1 = get(key1)
            def v2 = get(key2)
            def ret = v1 <=> v2
            if( ret == 0 ) {
                int i1 = ByteBuffer.wrap(key1).getInt()
                int i2 = ByteBuffer.wrap(key2).getInt()
                ret = i1 <=> i2
            }
            return ret
        }

        @Override
        public byte[] findShortestSeparator(byte[] start, byte[] limit) {
            return start;
        }

        @Override
        public byte[] findShortSuccessor(byte[] key) {
            return key;
        }

    }

    static class Cache<V> extends TCustomHashMap<byte[],V> {

        final int max

        Cache(int size) {
            super(new CustomArrayHashing(), size)
            this.max = size
        }

        public V put(byte[] key, V value) {
            if( size()==max && !containsKey(key) ) {
                removeAt(0)
            }
            super.put(key,value)
        }

    }

    static class CustomArrayHashing implements HashingStrategy {

        @Override
        int computeHashCode(Object o) {
            byte[] c = (byte[]) o;
            Arrays.hashCode(c)
        }

        @Override
        boolean equals(Object o1, Object o2) {
            byte[] c1 = (byte[])o1;
            byte[] c2 = (byte[])o2;

            if (c1.length != c2.length) {
                return false;
            }

            int i=0
            int len = c1.length
            for (; i < len; i++) {
                if (c1[i] != c2[i]) {
                    return false;
                }
            }
            return true;
        }
    }

    def TestLevelDb4() {

        cache = new Cache(1000)

        Options options = new Options()
                .createIfMissing(true)
        data = factory.open(new File("/Users/pditommaso/Downloads/db/data"), options);

        options = new Options()
                .createIfMissing(true)
                .comparator(new MyComparator())
        index = factory.open(new File("/Users/pditommaso/Downloads/db/index"), options);

    }

    def add( String line ) {
        //println "Adding $count"
        def key = bytes(count)
        cache.put(key, line)
        data.put(key, KryoHelper.serialize(line) )
        index.put(key, EMPTY)
        count++

        if( count % 5000 == 0 )
            println stats()
    }

    def stats() {
        "Added $count lines -- cache hits: $cacheHits - missed: $cacheMissed "
    }

    def close() {
        data.close()
        index.close()
    }

    def foreach(Closure<Void> closure) {
        DBIterator iterator = index.iterator();
        try {
            for(iterator.seekToFirst(); iterator.hasNext(); iterator.next()) {
                def key = (byte[])iterator.peekNext().getKey()
                def val = data.get(key)
                def ret = (String)KryoHelper.deserialize(val)
                closure(ret)
            }
        }
        finally {
            // Make sure you close the iterator to avoid resource leaks.
            iterator.close();
        }
    }

    void main(String[] args) {
        def file = new File('/Users/pditommaso/Downloads/db/')
        file.deleteDir()

        def collector = new TestLevelDb4()

        println "Begin"
        def now = System.currentTimeMillis()
        def text = new File('/Users/pditommaso/Downloads/hs_ref_GRCh38_chr1.fa') // hs_ref_GRCh38_chr1.fa
        text.eachLine { String it ->
            collector.add(it)
        }
        println "End adding - time: ${(System.currentTimeMillis()-now)/1000} secs -- Stats: ${collector.stats()}"


        println "\nSorted result: \n"
        def sort = new FileWriter(new File('/Users/pditommaso/Downloads/sort'))
        collector.foreach { sort.println it }
        sort.close()

        println "End - time: ${(System.currentTimeMillis()-now)/1000} secs"
        collector.close();

    }

}

