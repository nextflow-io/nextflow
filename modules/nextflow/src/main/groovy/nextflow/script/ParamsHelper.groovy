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

package nextflow.script

import java.lang.reflect.ParameterizedType
import java.lang.reflect.Type
import java.nio.file.Path
import java.util.function.BiFunction

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.yaml.YamlSlurper
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.dataflow.ChannelImpl
import nextflow.dataflow.ChannelNamespace
import nextflow.dataflow.ValueImpl
import nextflow.exception.ScriptRuntimeException
import nextflow.script.dsl.Nullable
import nextflow.script.dsl.PipelineParams
import nextflow.script.dsl.Types
import nextflow.script.types.Channel
import nextflow.script.types.Record
import nextflow.script.types.Value
import nextflow.splitter.CsvSplitter
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import nextflow.util.RecordMap
import nextflow.util.TypeHelper
import nextflow.util.VersionNumber
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
/**
 * Resolves a declared param against a given value.
 *
 * Used by pipeline params, typed workflows, and typed processes, so
 * that they all map a given value to a declared type the same way.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class ParamsHelper {

    /**
     * Resolve declared params from the command line and config.
     *
     * The config params already include the command line overrides
     * (see ConfigDsl), so they take precedence when present.
     *
     * @param declarations
     * @param cliParams
     * @param configParams
     */
    static Map<String,Object> resolveParams(Collection<Param> declarations, Map cliParams, Map configParams) {
        final names = declarations*.name as Set<String>
        for( final name : cliParams.keySet() ) {
            if( name !in names && !configParams.containsKey(name) )
                throw new ScriptRuntimeException("Parameter `${name}` was specified on the command line or params file but is not declared in the script or config")
        }

        final given = givenParams(names, cliParams, configParams)
        return resolveParams(declarations, given, '') { Param decl, Object value ->
            resolveParam(decl, value, cliParams.containsKey(decl.name))
        }
    }

    /**
     * Resolve declared params from the command line and config to plain
     * values, which can be serialized before the dataflow network has
     * started (e.g. for lineage or Seqera Platform).
     *
     * Each param is resolved as in {@link #resolveParams(Collection,Map,Map)},
     * except that a {@code Channel<E>} param is resolved to its samplesheet
     * and a {@code Value<V>} param to its value of type {@code V}, instead
     * of a dataflow value. The params are assumed to be valid, i.e. already
     * resolved by {@link #resolveParams(Collection,Map,Map)}.
     *
     * @param declarations
     * @param cliParams
     * @param configParams
     */
    static Map<String,Object> resolvePlainParams(Collection<Param> declarations, Map cliParams, Map configParams) {
        final given = givenParams(declarations*.name as Set<String>, cliParams, configParams)
        final result = new LinkedHashMap<String,Object>(declarations.size())
        for( final decl : declarations ) {
            final name = decl.name
            final value = given.containsKey(name)
                ? resolveParam(decl, given.get(name), cliParams.containsKey(name), true)
                : resolveDefault(decl, true)
            result.put(name, value)
        }
        return result
    }

    private static Map<String,?> givenParams(Set<String> names, Map cliParams, Map configParams) {
        return cliParams.subMap(names) + configParams.subMap(names)
    }

    /**
     * Resolve declared params against the given values. A param
     * with no given value is given its default value.
     *
     * @param declarations
     * @param given
     * @param context appended to the param name in error messages
     * @param resolve resolves a given value against its declared param
     */
    static Map<String,Object> resolveParams(Collection<Param> declarations, Map<String,?> given, String context, BiFunction<Param,Object,Object> resolve) {
        final result = new LinkedHashMap<String,Object>(declarations.size())
        for( final decl : declarations ) {
            final name = decl.name
            final value = given.containsKey(name)
                ? resolve.apply(decl, given.get(name))
                : resolveDefault(decl)

            if( value == null && !decl.optional )
                throw new ScriptRuntimeException("Parameter `${name}`${context} is required but no value was provided")

            result.put(name, value)
        }
        return result
    }

    /**
     * Resolve the params given to the entry workflow of a pipeline against
     * the params block of the pipeline. Called by the entry workflow (see
     * WorkflowToGroovyVisitor).
     *
     * The session params of a top-level run are already resolved by the
     * params block, so they are returned as-is. An included pipeline
     * resolves its params even when given the session params (e.g.
     * `RNASEQ(params)`), so that its defaults are applied.
     *
     * @param script the pipeline script
     * @param value the params given to the entry workflow
     */
    static Map resolveArguments(BaseScript script, Object value) {
        if( !ScriptMeta.get(script).isModule() )
            return (Map)value

        final pipeline = ExecutionStack.workflow().name
        if( value !instanceof RecordMap && value !instanceof ScriptBinding.ParamsMap )
            throw new ScriptRuntimeException("Pipeline `${pipeline}` should be called with a record")

        final given = (Map<String,?>)value
        final declarations = script.getParamDeclarations()
        for( final name : given.keySet() ) {
            if( !declarations.containsKey(name) )
                throw new ScriptRuntimeException("Pipeline `${pipeline}` does not declare a parameter named `${name}`")
        }

        final params = resolveParams(declarations.values(), given, " of pipeline `${pipeline}`") { Param decl, Object v ->
            isDataflow(v) ? DataflowTypeHelper.normalizeV2(v) : resolveParam(decl, v, false)
        }
        return new RecordMap(params)
    }

    private static boolean isDataflow(Object value) {
        return value instanceof ChannelImpl
            || value instanceof ValueImpl
            || value instanceof DataflowWriteChannel
            || value instanceof ChannelOut
    }

    /**
     * Resolve a param value against its declared type.
     *
     * A {@code Channel<E>} param is loaded from a samplesheet file, with each
     * record converted to the element type. A {@code Value<V>} param is
     * converted to {@code V} and wrapped in a dataflow value. Any other param
     * is converted directly to the declared type.
     *
     * @param decl
     * @param value
     * @param fromCli whether the value came from the command line (and is
     *                therefore a string that may need to be parsed)
     * @param plain whether to give the plain value of a {@code Channel<E>}
     *              or {@code Value<V>} param instead of a dataflow value
     *              (see {@link #resolvePlainParams})
     */
    static Object resolveParam(Param decl, Object value, boolean fromCli, boolean plain=false) {
        if( value == null )
            return null

        final rawType = TypeHelper.getRawType(decl.type)

        if( rawType == Channel )
            return plain ? value : ChannelNamespace.fromList(loadChannelInput(decl, value))

        if( rawType == Value ) {
            final result = resolveParam(elementDecl(decl), value, fromCli, plain)
            return plain ? result : ChannelNamespace.value(result)
        }

        if( TypeHelper.isRecordType(decl.type) && value instanceof Map )
            return resolveRecord(decl, (Map)value, fromCli, plain)

        final result = fromCli
            ? resolveFromCli(decl, value)
            : resolveFromCode(decl, value)
        checkAssignable(decl, result)
        return result
    }

    private static RecordMap resolveRecord(Param decl, Map value, boolean fromCli, boolean plain) {
        final type = (Class)decl.type
        final result = new LinkedHashMap<String,Object>(value)
        for( final field : type.getDeclaredFields() ) {
            if( field.isSynthetic() )
                continue
            final name = field.getName()
            final optional = field.isAnnotationPresent(Nullable)
            final fieldValue = value.get(name)
            if( fieldValue == null ) {
                if( !optional )
                    throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${type.getSimpleName()} is missing required field `${name}`")
                continue
            }
            final fieldDecl = new Param("${decl.name}.${name}", field.getGenericType(), optional, null)
            result.put(name, resolveParam(fieldDecl, fieldValue, fromCli, plain))
        }
        return new RecordMap(result)
    }

    /**
     * Load a channel param from a samplesheet file, converting each record
     * to the declared element type.
     *
     * @param decl
     * @param value
     */
    private static List loadChannelInput(Param decl, Object value) {
        if( value !instanceof CharSequence && value !instanceof Path )
            throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${Types.getName(decl.type)} should be a samplesheet file, but received: ${value} [${Types.getName(value.getClass())}]")

        final path = value instanceof Path
            ? (Path)value
            : TypeHelper.asPathType(value.toString())
        final elementType = elementDecl(decl).type
        final elementRawType = TypeHelper.getRawType(elementType)

        if( !Map.isAssignableFrom(elementRawType) && !Record.isAssignableFrom(elementRawType) )
            throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${Types.getName(decl.type)} cannot be loaded from a samplesheet -- the element type should be Map, Record, or a record type")

        return loadFromFile(decl.name, path).collect { el ->
            try {
                TypeHelper.asType(el, elementType)
            }
            catch( Exception e ) {
                throw new ScriptRuntimeException("Invalid record in samplesheet '${path}' for parameter `${decl.name}` -- ${e.message}")
            }
        }
    }

    /**
     * Get the declared param for the element type of a parameterized
     * type, e.g. {@code Sample} for {@code Channel<Sample>}.
     *
     * @param decl
     */
    private static Param elementDecl(Param decl) {
        final elementType = decl.type instanceof ParameterizedType
            ? ((ParameterizedType)decl.type).getActualTypeArguments()[0]
            : (Type)Object
        return new Param(decl.name, elementType, decl.optional, null)
    }

    /**
     * Load the contents of a samplesheet file as a list of records.
     *
     * Supported formats:
     * - CSV: header row required, comma-separated
     * - JSON: must be a top-level array
     * - YAML / YML: must be a top-level sequence
     *
     * @param name the param name (for error messages)
     * @param file the samplesheet file to load
     */
    static List loadFromFile(String name, Path file) {
        final ext = file.getExtension()
        final value = switch( ext ) {
            case 'csv'         -> loadFromCsv(file)
            case 'json'        -> new JsonSlurper().parse(file)
            case 'yaml', 'yml' -> new YamlSlurper().parse(file)
            default -> throw new ScriptRuntimeException("Unrecognized file format '${ext}' for input file '${file}' for parameter `${name}` -- should be CSV, JSON, or YAML")
        }
        if( value !instanceof List )
            throw new ScriptRuntimeException("Input file '${file}' for parameter `${name}` must contain a list of records, but got: ${value.class.simpleName}")
        return (List)value
    }

    private static List loadFromCsv(Path file) {
        final rows = new CsvSplitter().options(header: true, sep: ',', quote: '"').target(file).list()
        return rows.collect { row ->
            ((Map)row).collectEntries { k, v -> [ k, v != '' ? v : null ] }
        }
    }

    /**
     * Resolve a value given on the command line. Command-line values are
     * always strings, so they are parsed according to the declared type.
     *
     * @param decl
     * @param value
     */
    static Object resolveFromCli(Param decl, Object value) {
        if( value == null )
            return null

        if( value instanceof Collection || value instanceof Map )
            return asType(value, decl)

        final number = asNumberType(decl, value)
        if( number != null )
            return number

        if( value !instanceof CharSequence )
            return value

        final str = value.toString()

        if( decl.type == Boolean ) {
            if( str.toLowerCase() == 'true' ) return Boolean.TRUE
            if( str.toLowerCase() == 'false' ) return Boolean.FALSE
        }

        return resolveFromString(decl, str, value)
    }

    /**
     * Resolve a value given in a params file or the config. Such values are
     * already structured, so they only need to be converted where the
     * declared type is more specific than the source syntax.
     *
     * @param decl
     * @param value
     */
    static Object resolveFromCode(Param decl, Object value) {
        if( value == null )
            return null

        if( value instanceof Collection || value instanceof Map )
            return asType(value, decl)

        final number = asNumberType(decl, value)
        if( number != null )
            return number

        if( value !instanceof CharSequence )
            return value

        return resolveFromString(decl, value.toString(), value)
    }

    /**
     * Convert a value to a declared numeric type. Integer and Float are the
     * only numeric types in the Nextflow type system, but a value can be
     * given as any number (e.g. a Double or BigDecimal in the config), and a
     * command-line value is always a string, so the value is normalized to
     * the declared type.
     *
     * The value is converted to the narrowest representation that can hold
     * it -- Integer, Long, or BigInteger for an integral value, Float,
     * Double, or BigDecimal for a fractional one -- so that no precision is
     * lost for a value that is too large for the declared type.
     *
     * Returns null if the declared type is not numeric, or if the value
     * cannot be represented by it -- an Integer accepts an integral value in
     * any notation (e.g. `3`, `3.0`, `3e2`) but rejects one with a fractional
     * part, rather than silently truncating it. The value is then reported by
     * the caller as not assignable to the declared type.
     *
     * @param decl
     * @param value
     */
    private static Number asNumberType(Param decl, Object value) {
        if( decl.type != Integer && decl.type != Float )
            return null
        try {
            final number = new BigDecimal(value.toString().trim())
            return decl.type == Float
                ? asFloatType(number)
                : asIntegerType(number)
        }
        catch( NumberFormatException | ArithmeticException e ) {
            return null
        }
    }

    /**
     * Convert a number to a Float, widening it to a Double or BigDecimal
     * if it is too large to be represented by a Float.
     *
     * @param number
     */
    private static Number asFloatType(BigDecimal number) {
        final floatValue = number.floatValue()
        if( !floatValue.isInfinite() )
            return floatValue
        final doubleValue = number.doubleValue()
        return !doubleValue.isInfinite()
            ? doubleValue
            : number
    }

    /**
     * Convert a number to an Integer, widening it to a Long or BigInteger if
     * it is too large to be represented by an Integer. A fractional value is
     * rejected rather than being truncated.
     *
     * @param number
     */
    private static Number asIntegerType(BigDecimal number) {
        final integer = number.toBigIntegerExact()
        if( integer.bitLength() < Integer.SIZE )
            return integer.intValue()
        return integer.bitLength() < Long.SIZE
            ? integer.longValue()
            : integer
    }

    /**
     * Convert a string to a declared type that is always expressed as a
     * string, regardless of where the value came from. Returns the given
     * fallback if the declared type is not one of these.
     *
     * @param decl
     * @param str
     * @param fallback
     */
    private static Object resolveFromString(Param decl, String str, Object fallback) {
        if( decl.type == Path )
            return TypeHelper.asPathType(str)

        if( decl.type == Duration )
            return Duration.of(str)

        if( decl.type == MemoryUnit )
            return MemoryUnit.of(str)

        if( decl.type == VersionNumber )
            return new VersionNumber(str)

        return fallback
    }

    /**
     * Convert a composite value (a collection, map, or record) to the
     * declared type, reporting a conversion failure in terms of the param as
     * well as the underlying element -- a NumberFormatException from an
     * element of a `List<Integer>`, for example, names neither the param nor
     * the declared type on its own.
     *
     * @param value
     * @param decl
     */
    private static Object asType(Object value, Param decl) {
        try {
            return TypeHelper.asType(value, decl.type)
        }
        catch( GroovyCastException | UnsupportedOperationException | IllegalArgumentException e ) {
            final actualType = value.getClass()
            final detail = e.message ? " -- ${e.message}" : ''
            throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${value} [${Types.getName(actualType)}]${detail}")
        }
    }

    /**
     * The value of a param for which no value was provided: its
     * default value, if any, otherwise an empty record or null.
     *
     * A param whose type is the params block of an included pipeline
     * defaults to an empty record, so that the calling pipeline can
     * supply the params by dataflow. Missing params are reported when
     * the pipeline is called.
     *
     * @param decl
     * @param plain see {@link #resolveParam}
     */
    static Object resolveDefault(Param decl, boolean plain=false) {
        if( decl.defaultValue != null )
            return resolveParam(decl, decl.defaultValue, false, plain)
        final type = TypeHelper.getRawType(decl.type)
        return type.isAnnotationPresent(PipelineParams)
            ? new RecordMap([:])
            : null
    }

    /**
     * Check that a resolved value can be assigned to the declared type
     * of a param.
     *
     * @param decl
     * @param value
     */
    private static void checkAssignable(Param decl, Object value) {
        final expectedType = TypeHelper.getRawType(decl.type)
        final actualType = value?.getClass()
        if( actualType != null && !isAssignableFrom(expectedType, actualType) )
            throw new ScriptRuntimeException("Parameter `${decl.name}` with type ${Types.getName(decl.type)} cannot be assigned to ${value} [${Types.getName(actualType)}]")
    }

    static boolean isAssignableFrom(Class target, Class source) {
        // any numeric value can be assigned to Float
        if( target == Float.class )
            return Number.class.isAssignableFrom(source)

        // any integer value can be assigned to Integer
        if( target == Integer.class )
            return source == BigInteger.class || source == Long.class || source == Integer.class

        // any record can be assigned to a record type (validation is handled by asType())
        if( Record.class.isAssignableFrom(target) )
            return source == RecordMap.class

        return target.isAssignableFrom(source)
    }

}
