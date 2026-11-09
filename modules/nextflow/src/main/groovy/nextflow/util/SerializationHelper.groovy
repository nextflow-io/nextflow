/*
 * Copyright 2013-2026, Seqera Labs
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

import java.lang.reflect.Constructor
import java.lang.reflect.Modifier
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

import com.esotericsoftware.kryo.Kryo
import com.esotericsoftware.kryo.KryoException
import com.esotericsoftware.kryo.Serializer
import com.esotericsoftware.kryo.factories.SerializerFactory
import com.esotericsoftware.kryo.io.Input
import com.esotericsoftware.kryo.io.Output
import com.esotericsoftware.kryo.serializers.FieldSerializer
import de.javakaffee.kryoserializers.UnmodifiableCollectionsSerializer
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.io.SerializableMarker
import nextflow.plugin.Plugins
import org.codehaus.groovy.runtime.GStringImpl
import org.objenesis.instantiator.ObjectInstantiator
import org.objenesis.strategy.InstantiatorStrategy
import org.pf4j.Extension
/**
 * Helper class to get a {@code Kryo} object ready to be used
 */
@Slf4j
class KryoHelper {

    static final Map<Class,Object> serializers

    static final ThreadLocal<Kryo> threadLocal

    static {
        serializers = [:]

        for( SerializerRegistrant s : Plugins.getExtensions(SerializerRegistrant).iterator() ) {
            log.trace "Registering serializer: ${s.class.name}"
            s.register(serializers)
        }

        threadLocal = new ThreadLocal<Kryo>() {
            @Override
            protected Kryo initialValue() {
                newInstance()
            }
        }
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
        // use reflection/ASM instead of sun.misc.Unsafe in FieldSerializer
        kryo.setAsmEnabled(true)
        kryo.setInstantiatorStrategy( InstantiationStrategy.instance )
        kryo.setDefaultSerializer( DefaultSerializerFactory.instance )

        // special serializers
        // NOTE: classes registered here and below are written by ID, which
        // depends on the registration order -- changing it breaks existing caches
        UnmodifiableCollectionsSerializer.registerSerializers(kryo)

        // users serializers
        for( Map.Entry entry : serializers ) {
            final k = (Class)entry.key
            final v = entry.value
            if( v instanceof Class )
                kryo.register(k,(Serializer)v.newInstance())

            else if( v instanceof Closure )
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
        kryo.addDefaultSerializer(URI, URISerializer)

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
                kryo.setClassLoader(prev)
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

@Extension
@CompileStatic
class DefaultSerializers implements SerializerRegistrant {

    static private Class<Path> PATH_CLASS = (Class<Path>)Paths.get('.').class

    @Override
    void register(Map<Class, Object> serializers) {
        serializers.put( PATH_CLASS, PathSerializer )
        serializers.put( URL, URLSerializer )
        serializers.put( UUID, UUIDSerializer )
        serializers.put( File, FileSerializer )
        serializers.put( Pattern, PatternSerializer )
        serializers.put( ArrayTuple, ArrayTupleSerializer )
        serializers.put( RecordMap, RecordMapSerializer )
        serializers.put( SerializableMarker, null )
    }
}

/**
 * Creates objects through their no-arg constructor, without deep reflection into the JDK.
 *
 * A JDK collection or map without an accessible no-arg constructor (e.g. {@code List.of()},
 * {@code Arrays.asList()}) is created as the standard implementation of its interface.
 */
@Singleton
@CompileStatic
class InstantiationStrategy implements InstantiatorStrategy {

    @Override
    ObjectInstantiator newInstantiatorOf(final Class type) {
        final ctor = noArgConstructor(type)
        if( ctor != null )
            return new ObjectInstantiator() {
                @Override
                Object newInstance() {
                    try {
                        return ctor.newInstance()
                    }
                    catch( Exception e ) {
                        throw new KryoException("Error constructing instance of class: ${type.name}", e)
                    }
                }
            }
        if( Set.isAssignableFrom(type) )
            return { -> new LinkedHashSet() } as ObjectInstantiator
        if( Collection.isAssignableFrom(type) )
            return { -> new ArrayList() } as ObjectInstantiator
        if( Map.isAssignableFrom(type) )
            return { -> new LinkedHashMap() } as ObjectInstantiator
        throw new KryoException("Class cannot be created (missing no-arg constructor): ${type.name}")
    }

    static private Constructor noArgConstructor(Class type) {
        if( type.isMemberClass() && !Modifier.isStatic(type.getModifiers()) )
            return null
        try {
            final result = type.getDeclaredConstructor()
            return result.trySetAccessible() ? result : null
        }
        catch( NoSuchMethodException e ) {
            return null
        }
    }
}

/**
 * Serializes the fields of classes outside the JDK. A JDK type needs an explicit
 * serializer, since its fields are not accessible by reflection.
 */
@Singleton
@CompileStatic
class DefaultSerializerFactory implements SerializerFactory {

    @Override
    Serializer makeSerializer(Kryo kryo, Class<?> type) {
        if( type.getModule().isNamed() )
            throw new KryoException("No serializer is defined for type: ${type.name} -- only standard types can be cached")
        return new FieldSerializer(kryo, type)
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
    Path read(Kryo kryo, Input input, Class<Path> type) {
        final scheme = input.readString()
        final path = input.readString()
        log.trace "Path de-serialization > scheme: $scheme; path: $path"

        if( "file".equalsIgnoreCase(scheme) ) {
            return FileSystems.getDefault().getPath(path)
        }

        // NOTE: PATHS FOR NON DEFAULT FILE SYSTEMS SHOULD BE CREATED
        // BY THE CORRESPONDING SPECIALISED SERIALIZER CLASSES
        //
        // THIS IS LEGACY CODE THAT EVENTUALLY SHOULD BE REMOVED
        //
        def uri = URI.create("$scheme://$path")
        def fs = FileHelper.getOrCreateFileSystemFor(uri)
        return fs.provider().getPath(uri)
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

/**
 * Write a {@link URI} as its string, in the same format previously
 * written by the {@link FieldSerializer} for its only non-transient field
 */
@CompileStatic
class URISerializer extends Serializer<URI> {

    @Override
    void write(Kryo kryo, Output output, URI uri) {
        kryo.writeObjectOrNull(output, uri.toString(), String)
    }

    @Override
    URI read(Kryo kryo, Input input, Class<URI> type) {
        return new URI(kryo.readObjectOrNull(input, String))
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

@CompileStatic
class RecordMapSerializer extends Serializer<RecordMap> {

    @Override
    void write(Kryo kryo, Output output, RecordMap record) {
        final len = record.size()
        output.writeInt(len)
        for( Map.Entry<String,Object> entry : record.entrySet() ) {
            output.writeString(entry.getKey())
            kryo.writeClassAndObject(output, entry.getValue())
        }
    }

    @Override
    RecordMap read(Kryo kryo, Input input, Class<RecordMap> type) {
        final len = input.readInt()
        final props = new HashMap<String,Object>(len)
        for( int i=0; i<len; i++ ) {
            final key = input.readString()
            final value = kryo.readClassAndObject(input)
            props.put(key, value)
        }
        return new RecordMap(props)
    }
}
