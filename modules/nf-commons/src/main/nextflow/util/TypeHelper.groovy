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

import java.lang.reflect.Field
import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type
import java.nio.file.Files
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.script.dsl.Nullable
import nextflow.script.types.Bag
import nextflow.script.types.Record
import nextflow.util.HashBag
import nextflow.util.RecordMap
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
import org.codehaus.groovy.runtime.typehandling.GroovyCastException

/**
 * Utility functions for working with types at runtime.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TypeHelper {

    /**
     * Returns the given type argument from the given value's superclass.
     *
     * This method assumes that the value's class extends a parameterized type.
     *
     * @param value
     * @param index
     *
     * @example
     * <pre>
     * class ExampleClass extends GenericBase&lt;String, Integer&gt; {}
     *
     * Type type = TypeHelper.getGenericType(new ExampleClass(), 1);
     * System.out.println(type); // Output: class java.lang.Integer
     * </pre>
     */
    static Type getGenericType(Object value, int index) {
        final params = (ParameterizedType) value.getClass().getGenericSuperclass()
        return params.getActualTypeArguments()[index]
    }

    /**
     * Get the concrete Java class for a type.
     *
     * @param type
     */
    static Class getRawType(Type type) {
        return \
            type instanceof Class ? type :
            type instanceof ParameterizedType ? type.getRawType() :
            Object
    }

    /**
     * Determine whether a type is a collection type (Bag, List, Set).
     *
     * @param type
     */
    static boolean isCollectionType(Type type) {
        return Collection.class.isAssignableFrom(getRawType(type))
    }

    /**
     * Determine whether a type is a record type.
     *
     * @param type
     */
    static boolean isRecordType(Type type) {
        return type instanceof Class && Record.class.isAssignableFrom(type)
    }

    /**
     * Convert a value to the given type.
     *
     * @param value
     * @param type
     */
    static Object asType(Object value, Type type) {
        if( value == null )
            return null

        if( isCollectionType(type) )
            return asCollectionType(value as Collection, type)

        if( isRecordType(type) )
            return asRecordType(value as Map, (Class) type)

        if( type == Path )
            return TypeHelper.asPathType(value.toString())

        return DefaultTypeTransformation.castToType(value, getRawType(type))
    }

    /**
     * Convert a collection to the given collection type.
     *
     * If the type specifies an element type, each element in the
     * collection is converted to that type.
     *
     * @param collection
     * @param type
     */
    static Collection asCollectionType(Collection collection, Type type) {
        if( type instanceof ParameterizedType ) {
            final elementType = type.getActualTypeArguments()[0]
            collection = collection.collect { el -> asType(el, elementType) }
        }

        return switch( getRawType(type) ) {
            case Bag.class -> new HashBag<>(collection)
            case List.class -> collection as List
            case Set.class -> collection as Set
            default -> collection
        }
    }

    /**
     * Convert a string representing a file path to a Path.
     * Report an error if the path does not exist.
     *
     * @param str
     */
    static Path asPathType(String str) {
        final result = FileHelper.asPath(str)
        if( !Files.exists(result) )
            throw new AbortOperationException("Input file '${str}' does not exist")
        return result
    }

    /**
     * Convert a map to a record, validating it against the given
     * record type.
     *
     * @param map
     * @param type
     */
    static Record asRecordType(Map<String,?> map, Class type) {
        final fields = recordFields(type)

        for( final field : fields.values() ) {
            if( field.isAnnotationPresent(Nullable.class) )
                continue
            if( map.get(field.getName()) == null )
                throw new AbortOperationException("Input record ${map} is missing field '${field.getName()}' required by record type '${type.getSimpleName()}'")
        }

        final result = new HashMap<String,Object>(map.size())
        for( final entry : map.entrySet() ) {
            final name = entry.key
            final value = fields.containsKey(name)
                ? asType(entry.value, fields[name].getGenericType())
                : entry.value
            result.put(name, value)
        }
        return new RecordMap(result)
    }

    @Memoized
    private static Map<String,Field> recordFields(Class type) {
        final fields = type.getDeclaredFields()
        final result = new HashMap<String,Field>(fields.size())
        for( final field : fields ) {
            if( field.isSynthetic() )
                continue
            result.put(field.getName(), field)
        }
        return result
    }

}
