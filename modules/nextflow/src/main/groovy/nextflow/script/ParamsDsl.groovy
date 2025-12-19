/*
 * Copyright 2013-2024, Seqera Labs
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

import nextflow.script.types.Bag
import nextflow.util.ArrayBag

import java.nio.file.Path

import groovy.json.JsonSlurper
import groovy.yaml.YamlSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.exception.ScriptRuntimeException
import nextflow.script.types.Types
import nextflow.splitter.CsvSplitter
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import org.codehaus.groovy.runtime.typehandling.DefaultTypeTransformation
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
/**
 * Implements the DSL for defining workflow params
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsDsl {

    private Map<String,Param> declarations = [:]

    void declare(String name, Class type, boolean optional, Object defaultValue = null) {
        declarations[name] = new Param(name, type, optional, defaultValue)
    }

    void apply(Session session) {
        final cliParams = session.cliParams ?: [:]
        final configParams = session.configParams ?: [:]

        for( final name : cliParams.keySet() ) {
            if( !declarations.containsKey(name) && !configParams.containsKey(name) )
                throw new ScriptRuntimeException("Parameter `$name` was specified on the command line or params file but is not declared in the script or config")
        }

        final params = new HashMap<String,?>()
        for( final name : declarations.keySet() ) {
            final decl = declarations[name]
            if( cliParams.containsKey(name) ) {
                params[name] = resolveFromCli(decl, cliParams[name])
            }
            else if( configParams.containsKey(name) ) {
                params[name] = resolveFromCode(decl, configParams[name])
            }
            else if( decl.defaultValue != null ) {
                params[name] = resolveFromCode(decl, decl.defaultValue)
            }
            else if( decl.optional ) {
                params[name] = null
            }
            else {
                throw new ScriptRuntimeException("Parameter `$name` is required but was not specified on the command line, params file, or config")
            }

            final actualType = params[name]?.getClass()
            if( actualType != null && !isAssignableFrom(decl.type, actualType) )
                throw new ScriptRuntimeException("Parameter `$name` with type ${Types.getName(decl.type)} cannot be assigned to ${params[name]} [${Types.getName(actualType)}]")
        }

        session.binding.setParams(params, true)
    }

    private Object resolveFromCli(Param decl, Object value) {
        if( value == null )
            return null

        if( value !instanceof CharSequence )
            return value

        final str = value.toString()

        if( decl.type == Boolean ) {
            if( str.toLowerCase() == 'true' ) return Boolean.TRUE
            if( str.toLowerCase() == 'false' ) return Boolean.FALSE
        }

        if( decl.type == Integer || decl.type == Float ) {
            if( str.isInteger() ) return str.toInteger()
            if( str.isLong() ) return str.toLong()
            if( str.isBigInteger() ) return str.toBigInteger()
        }

        if( decl.type == Float ) {
            if( str.isFloat() ) return str.toFloat()
            if( str.isDouble() ) return str.toDouble()
            if( str.isBigDecimal() ) return str.toBigDecimal()
        }

        if( decl.type == Duration ) {
            return Duration.of(str)
        }

        if( decl.type == MemoryUnit ) {
            return MemoryUnit.of(str)
        }

        if( Iterable.class.isAssignableFrom(decl.type) || Map.class.isAssignableFrom(decl.type)) {
            return resolveFromFile(decl.name, decl.type, FileHelper.asPath(str))
        }

        if( decl.type == Path ) {
            return FileHelper.asPath(str)
        }

        return value
    }

    private Object resolveFromCode(Param decl, Object value) {
        if( value == null )
            return null

        if( value !instanceof CharSequence )
            return value

        final str = value.toString()

        if( Iterable.class.isAssignableFrom(decl.type) || Map.class.isAssignableFrom(decl.type))
            return resolveFromFile(decl.name, decl.type, FileHelper.asPath(str))

        if( decl.type == Path )
            return FileHelper.asPath(str)

        return value
    }

    private Object resolveFromFile(String name, Class type, Path file) {
        final ext = file.getExtension()
        final value = switch( ext ) {
            case 'csv' -> new CsvSplitter().options(header: true, sep: ',').target(file).list()
            case 'json' -> new JsonSlurper().parse(file)
            case 'yaml' -> new YamlSlurper().parse(file)
            case 'yml' -> new YamlSlurper().parse(file)
            default -> throw new ScriptRuntimeException("Unrecognized file format '${ext}' for input file '${file}' supplied for parameter `${name}` -- should be CSV, JSON, or YAML")
        }
        try {
            return DefaultTypeTransformation.castToType(value, type)
        }
        catch( GroovyCastException e ) {
            if( Bag.class.isAssignableFrom(type) || value instanceof Collection )
                return new ArrayBag(value)
            final actualType = value.getClass()
            throw new ScriptRuntimeException("Parameter `${name}` with type ${Types.getName(type)} cannot be assigned to contents of '${file}' [${Types.getName(actualType)}]")
        }
    }

    private boolean isAssignableFrom(Class target, Class source) {
        if( target == Float.class )
            return Number.class.isAssignableFrom(source)

        if( target == Integer.class )
            return source == BigInteger.class || source == Long.class || source == Integer.class

        return target.isAssignableFrom(source)
    }

    @Canonical
    private static class Param {
        String name
        Class type
        boolean optional
        Object defaultValue
    }

}
