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
package nextflow.agent

import java.lang.reflect.Field
import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type
import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.script.dsl.Nullable
import nextflow.script.types.Record

/**
 * Derives a portable JSON-schema {@link Map} by reflecting on an agent output
 * record {@link Class}. The resulting map mirrors the subset of JSON Schema used
 * as the LLM structured-output contract (see the nf-agent plugin which maps it
 * onto langchain4j's {@code JsonSchema}).
 *
 * A named record type used as an agent output lowers to a concrete runtime class
 * whose {@code getDeclaredFields()} returns the declared fields; optional fields
 * (declared with the {@code ?} suffix) carry a {@link Nullable} annotation and are
 * omitted from the {@code required} list.
 *
 * Supported field types (v1): String, integer/long, floating point/decimal,
 * boolean, nested record types, and List/Collection/Set of those. {@link Path}
 * and any other unmapped type are rejected with an {@link IllegalArgumentException}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class RecordSchema {

    /**
     * Reflect the given record type into a JSON-schema map.
     *
     * @param recordType the compiled record class (e.g. an agent output type)
     * @return an ordered map with {@code type}, {@code properties}, {@code required}
     *         and {@code additionalProperties} entries
     */
    static Map of(Class recordType) {
        final properties = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()

        for( final field : recordType.getDeclaredFields() ) {
            if( field.isSynthetic() )
                continue
            final name = field.getName()
            properties.put(name, fragmentFor(name, field.getGenericType()))
            if( !field.isAnnotationPresent(Nullable.class) )
                required.add(name)
        }

        final result = new LinkedHashMap<String,Object>()
        result.put('type', 'object')
        result.put('properties', properties)
        result.put('required', required)
        result.put('additionalProperties', false)
        return result
    }

    private static Map fragmentFor(String fieldName, Type type) {
        final raw = rawClass(type)

        if( raw == String )
            return [type: 'string']

        if( raw in [Integer, int, Long, long, Short, short, Byte, byte, BigInteger] )
            return [type: 'integer']

        if( raw in [Double, double, Float, float, BigDecimal, Number] )
            return [type: 'number']

        if( raw in [Boolean, boolean] )
            return [type: 'boolean']

        if( raw == Path || Path.isAssignableFrom(raw) )
            throw new IllegalArgumentException("Unsupported agent output field `${fieldName}` of type ${raw.getName()} - `Path` is not allowed in agent outputs")

        if( Record.isAssignableFrom(raw) )
            return of(raw)

        if( Collection.isAssignableFrom(raw) ) {
            final elementType = elementType(type)
            return [type: 'array', items: fragmentFor("${fieldName}[]".toString(), elementType)]
        }

        throw new IllegalArgumentException("Unsupported agent output field `${fieldName}` of type ${raw.getName()} - supported types are String, integer, number, boolean, nested record and list of those")
    }

    private static Class rawClass(Type type) {
        if( type instanceof Class )
            return type
        if( type instanceof ParameterizedType )
            return (Class) type.getRawType()
        return Object
    }

    private static Type elementType(Type type) {
        if( type instanceof ParameterizedType ) {
            final args = type.getActualTypeArguments()
            if( args.length > 0 )
                return args[0]
        }
        return String
    }

}
