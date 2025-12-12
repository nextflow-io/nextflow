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

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.file.FileHelper
import nextflow.exception.ScriptRuntimeException
import nextflow.script.types.Types
/**
 * Implements the DSL for defining workflow params
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsDsl {

    private Map<String,Param> declarations = [:]

    void declare(String name, Class type) {
        declarations[name] = new Param(name, type, Optional.empty())
    }

    void declare(String name, Class type, Object defaultValue) {
        declarations[name] = new Param(name, type, Optional.of(defaultValue))
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
            if( cliParams.containsKey(name) )
                params[name] = resolveFromCli(decl, cliParams[name])
            else if( configParams.containsKey(name) )
                params[name] = resolveFromCode(decl, configParams[name])
            else if( decl.defaultValue.isPresent() )
                params[name] = resolveFromCode(decl, decl.defaultValue.get())
            else
                throw new ScriptRuntimeException("Parameter `$name` is required but was not specified on the command line, params file, or config")

            final actualType = params[name].getClass()
            if( !isAssignableFrom(decl.type, actualType) )
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

        if( decl.type == Path ) {
            return FileHelper.asPath(str)
        }

        return value
    }

    private Object resolveFromCode(Param decl, Object value) {
        if( value == null )
            return null

        if( decl.type == Path && value instanceof CharSequence )
            return FileHelper.asPath(value.toString())

        return value
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
        Optional<?> defaultValue
    }

}
