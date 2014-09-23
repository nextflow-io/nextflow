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

import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import groovy.util.logging.Slf4j
import nextflow.io.DataInputStreamAdapter
import nextflow.io.DataOutputStreamAdapter
import nextflow.util.KryoHelper
import org.mapdb.BTreeKeySerializer
import org.mapdb.DB
import org.mapdb.DBMaker
import org.mapdb.Serializer
/**
 *  Helper class used to aggregate values having the same key
 *  to files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FileCollector implements Closeable {

    /**
     * Hold a single entry in the tree-set
     */
    static class Entry implements Serializable {

        Comparable name
        Comparable value
        long index

        Entry( Comparable a, Comparable b, long n ) {
            this.name = a
            this.value = b
            this.index = n
        }

    }

    /**
     * Compare entries using the {@link Entry#name} as primary index and
     * the value itself and secondary key
     */
    static class EntryComp implements Comparator<Entry>, Serializable {

        transient Closure<Comparable> sort

        @Override
        int compare(Entry e1, Entry e2) {

            def k = e1.name <=> e2.name
            if( k != 0 )
                return k

            if( sort )
                return sort.call(e1.value) <=> sort.call(e2.value)

            return e1.index <=> e2.index
        }
    }

    /**
     * Serialize an entry as required by MapDB
     */
    static class EntrySerializer implements Serializer<Entry>, Serializable {

        @Override
        void serialize(DataOutput out, Entry value) throws IOException {
            def output = new Output(new DataOutputStreamAdapter(out))
            KryoHelper.kryo().writeObject(output, value)
            output.flush()
        }

        @Override
        Entry deserialize(DataInput dataInput, int available) throws IOException {
            def input = new Input(new DataInputStreamAdapter(dataInput))
            return KryoHelper.kryo().readObject(input,Entry)
        }

        @Override
        int fixedSize() { -1 }
    }

    static final OpenOption[] APPEND = [StandardOpenOption.CREATE, StandardOpenOption.APPEND, StandardOpenOption.WRITE] as OpenOption[]

    public Boolean newLine

    public seed

    public Closure sort

    private DB db

    private SortedSet<Entry> cache

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


    private SortedSet<Entry> getOrCreateCache() {

        if( cache )
            return cache

        this.db = DBMaker.newFileDB(tempDbFile)
                .transactionDisable()
                .asyncWriteEnable()
                .mmapFileEnable()
                .compressionEnable()
                .closeOnJvmShutdown()
                .deleteFilesAfterClose()
                .make()

        def rndName = "set-${new BigInteger(130, new Random()).toString(32)}"
        cache = db
                .createTreeSet(rndName)
                .comparator(new EntryComp(sort: sort))
                .serializer(new BTreeKeySerializer.BasicKeySerializer(new EntrySerializer()))
                .makeOrGet()

        return cache
    }

    protected InputStream source( value ) {
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


    FileCollector append( String key, Comparable value ) {
        getOrCreateCache().add( new Entry(key,value,count++) )
        return this
    }


    /**
     * Append the content of a file to the target file having {@code key} as name
     *
     * @param key
     * @param fileToAppend
     */
    protected void append( InputStream source, OutputStream target ) {
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

        if( !cache ) {
            // nothing  to do
            return
        }

        def last = null
        OutputStream output = null

        Iterator<Entry> itr = cache.iterator()
        while( itr.hasNext() ) {
            def entry = itr.next()
            def name = entry.name
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
                    append(source(seed.get(name)), output)
                }
                else if( seed ) {
                    append(source(seed), output)
                }

            }

            /*
             * add the current value
             */
            append(source(entry.value), output)

        }

        // close last output stream
        output?.closeQuietly()
    }

    @Override
    void close() {
        db?.close()
    }


}
