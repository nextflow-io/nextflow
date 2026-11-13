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
package de.javakaffee.kryoserializers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;

/**
 * A kryo {@link Serializer} for unmodifiable {@link Collection}s and {@link Map}s
 * created via {@link Collections}.
 *
 * The wrapped collection is written as a copy, instead of reading the private
 * field of the wrapper, so that no reflective access to {@code java.util} is required.
 *
 * @author Martin Grotzke <martin.grotzke@javakaffee.de>
 *
 * @link https://github.com/magro/kryo-serializers
 */
public class UnmodifiableCollectionsSerializer extends Serializer<Object> {

    @Override
    public Object read(final Kryo kryo, final Input input, final Class<Object> clazz) {
        final int ordinal = input.readInt( true );
        final UnmodifiableCollection unmodifiableCollection = UnmodifiableCollection.values()[ordinal];
        final Object sourceCollection = kryo.readClassAndObject( input );
        return unmodifiableCollection.create( sourceCollection );
    }

    @Override
    public void write(final Kryo kryo, final Output output, final Object object) {
        final UnmodifiableCollection unmodifiableCollection = UnmodifiableCollection.valueOfType( object.getClass() );
        // the ordinal could be replaced by s.th. else (e.g. a explicitly managed "id")
        output.writeInt( unmodifiableCollection.ordinal(), true );
        kryo.writeClassAndObject( output, unmodifiableCollection.copy( object ) );
    }

    private static enum UnmodifiableCollection {
        COLLECTION( Collections.unmodifiableCollection( Arrays.asList( "" ) ).getClass() ){
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableCollection( (Collection<?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new ArrayList<>( (Collection<?>) collection );
            }
        },
        RANDOM_ACCESS_LIST( Collections.unmodifiableList( new ArrayList<Void>() ).getClass() ){
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableList( (List<?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new ArrayList<>( (Collection<?>) collection );
            }
        },
        LIST( Collections.unmodifiableList( new LinkedList<Void>() ).getClass() ){
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableList( (List<?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new ArrayList<>( (Collection<?>) collection );
            }
        },
        SET( Collections.unmodifiableSet( new HashSet<Void>() ).getClass() ){
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableSet( (Set<?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new LinkedHashSet<>( (Collection<?>) collection );
            }
        },
        SORTED_SET( Collections.unmodifiableSortedSet( new TreeSet<Void>() ).getClass() ){
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableSortedSet( (SortedSet<?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new TreeSet<>( (SortedSet<?>) collection );
            }
        },
        MAP( Collections.unmodifiableMap( new HashMap<Void, Void>() ).getClass() ) {
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableMap( (Map<?, ?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new LinkedHashMap<>( (Map<?, ?>) collection );
            }
        },
        SORTED_MAP( Collections.unmodifiableSortedMap( new TreeMap<Void, Void>() ).getClass() ) {
            @Override
            public Object create( final Object sourceCollection ) {
                return Collections.unmodifiableSortedMap( (SortedMap<?, ?>) sourceCollection );
            }
            @Override
            public Object copy( final Object collection ) {
                return new TreeMap<>( (SortedMap<?, ?>) collection );
            }
        };

        private final Class<?> type;

        private UnmodifiableCollection( final Class<?> type ) {
            this.type = type;
        }

        /**
         * @param sourceCollection
         */
        public abstract Object create( Object sourceCollection );

        /**
         * @param collection the unmodifiable collection
         * @return a modifiable copy of the given collection
         */
        public abstract Object copy( Object collection );

        static UnmodifiableCollection valueOfType( final Class<?> type ) {
            for( final UnmodifiableCollection item : values() ) {
                if ( item.type.equals( type ) ) {
                    return item;
                }
            }
            throw new IllegalArgumentException( "The type " + type + " is not supported." );
        }

    }

    /**
     * Creates a new {@link UnmodifiableCollectionsSerializer} and registers its serializer
     * for the several unmodifiable Collections that can be created via {@link Collections},
     * including {@link Map}s.
     *
     * NOTE: the registration order determines the class IDs written in the cache, so it must not change.
     *
     * @param kryo the {@link Kryo} instance to set the serializer on.
     *
     * @see Collections#unmodifiableCollection(Collection)
     * @see Collections#unmodifiableList(List)
     * @see Collections#unmodifiableSet(Set)
     * @see Collections#unmodifiableSortedSet(SortedSet)
     * @see Collections#unmodifiableMap(Map)
     * @see Collections#unmodifiableSortedMap(SortedMap)
     */
    public static void registerSerializers( final Kryo kryo ) {
        final UnmodifiableCollectionsSerializer serializer = new UnmodifiableCollectionsSerializer();
        for ( final UnmodifiableCollection item : UnmodifiableCollection.values() ) {
            kryo.register( item.type, serializer );
        }
    }

}
