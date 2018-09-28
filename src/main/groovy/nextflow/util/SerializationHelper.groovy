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

package nextflow.util
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Matcher
import java.util.regex.Pattern

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.upplication.s3fs.S3Path
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.file.http.XPath
import org.codehaus.groovy.runtime.GStringImpl
import org.objenesis.instantiator.ObjectInstantiator

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
        serializers.put( UUID, UUIDSerializer )
        serializers.put( File, FileSerializer )
        serializers.put( S3Path, PathSerializer )
        serializers.put( XPath, XPathSerializer )
        serializers.put( Pattern, PatternSerializer )
        serializers.put( ArrayTuple, ArrayTupleSerializer )

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
        kryo.setInstantiatorStrategy( InstantiationStrategy.instance )

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

        // Register default serializers
        // - there are multiple GString implementation so it
        //   has to be registered the base class as default serializer
        //   See Nextflow issues #12
        kryo.addDefaultSerializer(GString, GStringSerializer)
        // map entry serializer
        kryo.addDefaultSerializer(Map.Entry, MapEntrySerializer)

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
 * Kryo throws an IllegalAccessError when deserializing a {@link Matcher} object,
 * thus implements a custom instantiator object
 */
@Singleton
@CompileStatic
class MatcherInstantiator implements ObjectInstantiator {

    @Override
    Object newInstance() {
        def c = Matcher.getDeclaredConstructor()
        c.setAccessible(true)
        c.newInstance()
    }
}

@Singleton
@CompileStatic
class InstantiationStrategy extends Kryo.DefaultInstantiatorStrategy {

    @Override
    public ObjectInstantiator newInstantiatorOf (final Class type) {
        if( type == Matcher ) {
            MatcherInstantiator.instance
        }
        else {
            super.newInstantiatorOf(type)
        }
    }
}

/**
 * A Kryo serializer to handler a {@code Path} object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@CompileStatic
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

@Slf4j
@CompileStatic
class XPathSerializer extends Serializer<Path> {

    @Override
    void write(Kryo kryo, Output output, Path target) {
        final uri = target.toUri().toString()
        log.trace "XPath serialization > uri=$uri"
        output.writeString(uri)
    }

    @Override
    Path read(Kryo kryo, Input input, Class<Path> type) {
        final uri = input.readString()
        log.trace "Path de-serialization > uri=$uri"
        FileHelper.asPath(new URI(uri))
    }
}

/**
 * Serializer / de-serializer for Groovy GString
 */
@Slf4j
@CompileStatic
class GStringSerializer extends Serializer<GString> {

    static final private Class<Object[]> OBJ_ARRAY_CLASS = (Class<Object[]>)(new Object[0]).getClass()
    static final private Class<String[]> STR_ARRAY_CLASS = (Class<String[]>)(new String[0]).getClass()

    @Override
    void write(Kryo kryo, Output stream, GString object) {
        log.trace "GString serialization: values: ${object?.getValues()} - strings: ${object?.getStrings()}"
        kryo.writeObject( stream, object.getValues() )
        kryo.writeObject( stream, object.getStrings() )
    }

    @Override
    GString read(Kryo kryo, Input stream, Class<GString> type) {
        Object[] values = kryo.readObject(stream, OBJ_ARRAY_CLASS)
        String[] strings = kryo.readObject(stream, STR_ARRAY_CLASS)
        log.trace "GString de-serialize: values: ${values} - strings: ${strings}"
        new GStringImpl(values, strings)
    }
}


@Slf4j
@CompileStatic
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
@CompileStatic
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
@CompileStatic
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


/**
 * The {@link Pattern} class does not have default constructor as required by Kryo,
 * but implements the {@link Serializable} interface, thus use the default serialization
 * mechanism to store the pattern as a byte array
 */
@CompileStatic
class PatternSerializer extends Serializer<Pattern> {

    @Override
    void write(Kryo kryo, Output output, Pattern object) {

        def buffer = new ByteArrayOutputStream()
        ObjectOutputStream oos = new ObjectOutputStream(buffer)
        oos.writeObject(object)
        oos.close()

        def bytes = buffer.toByteArray()
        output.writeInt(bytes.length)
        output.write(bytes)
    }

    @Override
    Pattern read(Kryo kryo, Input input, Class<Pattern> type) {

        def len = input.readInt()
        def buffer = new byte[len]
        input.read(buffer)

        ObjectInputStream ois = new ObjectInputStream(new ByteArrayInputStream(buffer))
        return (Pattern) ois.readObject()

    }
}

@CompileStatic
class ArrayTupleSerializer extends Serializer<ArrayTuple> {

    @Override
    void write(Kryo kryo, Output output, ArrayTuple tuple) {
        final len = tuple.size()
        output.writeInt(len)
        for( int i=0; i<len; i++ ) {
            kryo.writeClassAndObject(output, tuple.get(i))
        }
    }

    @Override
    ArrayTuple read(Kryo kryo, Input input, Class<ArrayTuple> type) {
        final len = input.readInt()
        def list = new ArrayList(len)
        for( int i=0; i<len; i++ ) {
            def item = kryo.readClassAndObject(input)
            list.add(item)
        }
        return new ArrayTuple(list)
    }
}

@CompileStatic
class MapEntrySerializer extends Serializer<Map.Entry> {

    @Override
    void write(Kryo kryo, Output output, Map.Entry entry) {
        kryo.writeClassAndObject(output, entry.getKey())
        kryo.writeClassAndObject(output, entry.getValue())
    }

    @Override
    Map.Entry read(Kryo kryo, Input input, Class<Map.Entry> type) {
        def key = kryo.readClassAndObject(input)
        def val = kryo.readClassAndObject(input)
        new MapEntry(key,val)
    }
}
