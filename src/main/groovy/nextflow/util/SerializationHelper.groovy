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

package nextflow.util
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer
import groovy.util.logging.Slf4j
import org.codehaus.groovy.runtime.GStringImpl
/**
 * Helper class to get a {@code Kryo} object ready to be used
 */
class KryoHelper {

    static private Class<Path> PATH_CLASS = Paths.get('.').class

    static final Map<Class,Object> serializers

    static final ThreadLocal<Kryo> threadLocal

    static {
        serializers = [:]
        serializers.put( PATH_CLASS, PathSerializer )
        serializers.put( URL, URLSerializer )
        serializers.put( GStringImpl, GStringSerializer)
        serializers.put( UUID, UUIDSerializer )
        serializers.put( File, FileSerializer )

        threadLocal = new ThreadLocal<Kryo>() {
            @Override
            protected Kryo initialValue() {
                newInstance()
            }
        };
    }

    /**
     * Register a new class - serializer pair
     *
     * @param clazz
     * @param item
     */
    static void register( Class clazz, item = null ) {
        serializers.put( clazz, item )
    }

    /**
     * @return A new instance {@code Kryo} instance
     */
    static Kryo kryo() {
        threadLocal.get()
    }


    /**
     * @return A new instance {@code Kryo} instance
     */
    static private Kryo newInstance() {
        def kryo = new Kryo()
        // special serializers
        UnmodifiableCollectionsSerializer.registerSerializers(kryo)

        // users serializers
        serializers.each { k, v ->
            if( v instanceof Class )
                kryo.register(k,(Serializer)v.newInstance())

            else if( v instanceof Closure<Serializer> )
                kryo.register(k, v.call(kryo))

            else
                kryo.register(k)
        }

        return kryo
    }

    /**
     * De-serialize an object stored to a file
     *
     * @param file A {@code Path} reference to a file
     * @return The de-serialized object
     */
    static <T> T deserialize( Path file ) {
        Input input = null
        try {
            input = new Input(Files.newInputStream(file))
            return (T)kryo().readClassAndObject(input)
        }
        finally {
            input?.closeQuietly()
        }
    }

    /**
     * De-serialize an object stored to a file
     *
     * @param file A {@code File} reference to a file
     * @return The de-serialized object
     */
    static <T> T deserialize( File file ) {
        (T)deserialize(file.toPath())
    }

    /**
     * Serialize an object to the specified file
     *
     * @param object The object to serialize
     * @param toFile A {@code Path} reference to the file where store the object
     */
    static void serialize( object, Path toFile ) {
        def output = null
        try {
            output = new Output(Files.newOutputStream(toFile))
            kryo().writeClassAndObject(output, object)
        }
        finally {
            output?.closeQuietly()
        }
    }

    static byte[] serialize( object ) {
        ByteArrayOutputStream buffer = new ByteArrayOutputStream(4*1024)
        def output = new Output(buffer)
        kryo().writeClassAndObject(output, object)
        output.flush()
        return buffer.toByteArray()
    }

    static <T> T deserialize( byte[] binary, ClassLoader loader = null ) {
        def kryo = kryo()
        def ClassLoader prev = null
        if( loader ) {
            prev = kryo.getClassLoader()
            kryo.setClassLoader(loader)
        }

        try {
            def buffer = new ByteArrayInputStream(binary)
            return (T)kryo.readClassAndObject(new Input(buffer))
        }
        finally {
            if( prev ) {
                kryo.setClassLoader(loader)
            }
        }
    }


    /**
     * Serialize an object to the specified file
     *
     * @param object The object to serialize
     * @param toFile A {@code File} reference to the file where store the object
     */
    static void serialize( object, File toFile ) {
        serialize(object,toFile.toPath())
    }

}


/**
 * A Kryo serializer to handler a {@code Path} object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
class PathSerializer extends Serializer<Path> {

    @Override
    void write(Kryo kryo, Output output, Path target) {
        final scheme = target.getFileSystem().provider().getScheme()
        final path = target.toString()
        log.trace "Path serialization > scheme: $scheme; path: $path"
        output.writeString(scheme)
        output.writeString(path)
    }

    @Override
    Path  read(Kryo kryo, Input input, Class<Path> type) {
        final scheme = input.readString()
        final path = input.readString()
        log.trace "Path de-serialization > scheme: $scheme; path: $path"

        if( "file".equalsIgnoreCase(scheme) ) {
            return FileSystems.getDefault().getPath(path)
        }

        // try to find provider
        def uri = URI.create("$scheme://$path")
        def fs = FileHelper.getOrCreateFileSystemFor(uri)
        return fs.provider().getPath(uri)

    }
}

/**
 * Serializer / de-serializer for Groovy GString
 */
@Slf4j
class GStringSerializer extends Serializer<GStringImpl> {

    @Override
    void write(Kryo kryo, Output stream, GStringImpl object) {
        kryo.writeObject( stream, object.getValues() )
        kryo.writeObject( stream, object.getStrings() )
    }

    @Override
    GStringImpl read(Kryo kryo, Input stream, Class<GStringImpl> type) {
        Object[] values = kryo.readObject(stream, Object[].class)
        String[] strings = kryo.readObject(stream, String[].class)
        return new GStringImpl(values, strings)
    }
}



@Slf4j
class URLSerializer extends Serializer<URL> {

    @Override
    void write(Kryo kryo, Output output, URL url) {
        log.trace "URL serialization > $url"
        output.writeString(url.toString())
    }

    @Override
    URL read(Kryo kryo, Input input, Class<URL> type) {
        log.trace "URL de-serialization"
        return new URL(input.readString())
    }
}

@Slf4j
class UUIDSerializer extends Serializer<UUID> {

    @Override
    void write(Kryo kryo, Output output, UUID obj) {
        log.trace "UUID serialization > $obj"
        output.writeLong( obj.mostSignificantBits )
        output.writeLong( obj.leastSignificantBits )
    }

    @Override
    UUID read(Kryo kryo, Input input, Class<UUID> type) {
        log.trace "UUID de-serialization"
        long mostBits = input.readLong()
        long leastBits = input.readLong()
        return new UUID(mostBits, leastBits)
    }
}


@Slf4j
class FileSerializer extends Serializer<File> {

    @Override
    void write(Kryo kryo, Output output, File file) {
        log.trace "File serialization > $file"
        output.writeString(file.toString())
    }

    @Override
    File read(Kryo kryo, Input input, Class<File> type) {
        log.trace "File de-serialization"
        return new File(input.readString())
    }
}


