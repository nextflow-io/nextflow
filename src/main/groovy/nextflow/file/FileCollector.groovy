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
import java.nio.file.attribute.BasicFileAttributes

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
//import org.mapdb.DB
//import org.mapdb.DBMaker
//import org.mapdb.DataInput2
//import org.mapdb.Serializer

/**
 *  Helper class used to aggregate values having the same key
 *  to files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileCollector implements Closeable {

    /**
     * Hold a single entry in the tree-set
     */
    static class IndexEntry implements Serializable {

        /*
         * User provided grouping key used to sort entries in the (secondary) index
         */
        Comparable groupKey

        /*
         * The entry index i.e. its key in the store map
         */
        long index

        // required by kryo de-serialization
        protected IndexEntry () {}

        IndexEntry( Comparable a, long n ) {
            this.groupKey = a
            this.index = n
        }

    }

    /**
     * Compare entries using the {@link IndexEntry#groupKey} as primary index and
     * the value itself and secondary key
     */
    static class EntryComp implements Comparator<IndexEntry>, Serializable {

        transient Closure<Comparable> sort
        transient Map<Long,Comparable> store

        @Override
        int compare(IndexEntry e1, IndexEntry e2) {

            def k = e1.groupKey <=> e2.groupKey
            if( k != 0 )
                return k

            if( sort )
                return sort.call(store[e1.index]) <=> sort.call(store[e2.index])

            return e1.index <=> e2.index
        }
    }

    /**
     * Serialize an entry as required by MapDB
     */
//    static class EntrySerializer implements Serializer<IndexEntry>, Serializable {
//
//        @Override
//        void serialize(DataOutput dataOutput, IndexEntry value) throws IOException {
//            def output = new Output(dataOutput as OutputStream)
//            KryoHelper.kryo().writeObject(output, value)
//            output.flush()
//        }
//
//        @Override
//        IndexEntry deserialize(DataInput dataInput, int available) throws IOException {
//            def input = new Input(dataInput as InputStream)
//            return KryoHelper.kryo().readObject(input,IndexEntry)
//        }
//
//        @Override
//        int fixedSize() { -1 }
//    }

    /**
     * Serialize an entry as required by MapDB
     */
//    static class ObjectSerializer implements Serializer<Object>, Serializable {
//
//        @Override
//        void serialize(DataOutput dataOutput, Object value) throws IOException {
//            def output = new Output(dataOutput as OutputStream)
//            KryoHelper.kryo().writeClassAndObject(output, value)
//            output.flush()
//        }
//
//        @Override
//        Object deserialize(DataInput dataInput, int available) throws IOException {
//            def input = new Input(new DataInputPatch(dataInput as DataInput2))
//            return KryoHelper.kryo().readClassAndObject(input)
//        }
//
//        @Override
//        int fixedSize() { return -1 }
//    }

    static final OpenOption[] APPEND = [StandardOpenOption.CREATE, StandardOpenOption.APPEND, StandardOpenOption.WRITE] as OpenOption[]

    public Boolean newLine

    public seed

    public Closure sort

    //private DB db

    private Map<Long,Comparable> store

    private NavigableSet<IndexEntry> index

    private long count

    private File tempDbFile

    /**
     * Creates a file collector
     *
     * @param path The path the temporary files will be stored
     */
    FileCollector( Path path = null ) {
        if( !path ) {
            tempDbFile = Files.createTempFile('nxf', 'collect').toFile()
        }
        else {
            // read the path attributes
            def attr = Files.readAttributes(path, BasicFileAttributes)
            // when it is a dir => create a temp file there
            if ( attr.isDirectory() ) {
                tempDbFile = Files.createTempFile(path,'nxf','collect').toFile()
            }
            // if its a file use it
            else if( attr.isRegularFile()) {
                tempDbFile = path.toFile()
            }
            else {
                def dir = Files.createDirectories(path)
                tempDbFile = Files.createTempFile(dir,'nxf','collect').toFile()
            }

        }
    }

    /**
     * Creates the MapDB data structures
     */
    private void createStoreAndIndex() {

        /*
         * Creates the MapDB object
         */
//        db = DBMaker.newFileDB(tempDbFile)
//                .transactionDisable()
//                //.cacheDisable()
//                .mmapFileEnable()
//                .compressionEnable()
//                .closeOnJvmShutdown()
//                .deleteFilesAfterClose()
//                .make()

        def unique = new BigInteger(130, new Random()).toString(32)

        /*
         * The map that holds the objects
         */
//        store = db
//                .createHashMap("map-${unique}")
//                .keySerializer( Serializer.LONG )
//                .valueSerializer( new ObjectSerializer() )
//                .makeOrGet()

        /*
         * a "sorted set" used as index to access the values in
         */
//        index = db
//                .createTreeSet("ndx-${unique}")
//                .comparator(new EntryComp(sort: sort, store: store))
//                //.serializer(BTreeKeySerializer.TUPLE2)
//                .makeOrGet()

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
     * @return The {@link FileCollector} self object
     */
    FileCollector append( String key, Comparable value ) {
        if( !store )
            createStoreAndIndex()

        store.put( count, value )
        index.add( new IndexEntry(key,count) )
        count++

        return this
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

        Iterator<IndexEntry> itr = index.iterator()
        while( itr.hasNext() ) {
            def entry = itr.next()
            def name = entry.groupKey
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
                if( seed instanceof Map && seed.containsKey(name)) {
                    appendStream(normalizeToStream(seed.get(name)), output)
                }
                else if( seed ) {
                    appendStream(normalizeToStream(seed), output)
                }

            }

            /*
             * add the current value
             */
            def val = store.get(entry.index)
            appendStream(normalizeToStream(val), output)

        }

        // close last output stream
        output?.closeQuietly()
    }

    @Override
    void close() {
        //db?.close()
    }


}
