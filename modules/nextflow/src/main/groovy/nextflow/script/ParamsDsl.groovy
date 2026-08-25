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

import java.lang.reflect.Type

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.script.dsl.Types
import nextflow.util.TypeHelper
/**
 * Implements the DSL for defining workflow params
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsDsl {

    private Class clazz

    private Map<String,Param> declarations = [:]

    ParamsDsl(Class clazz) {
        this.clazz = clazz
    }

    void declare(String name, boolean optional, Object defaultValue = null) {
        final type = clazz.getField(name).getGenericType()
        if( defaultValue == null && type == Boolean )
            defaultValue = false
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
                params[name] = ParamsHelper.resolveFromCli(decl, cliParams[name])
            }
            else if( configParams.containsKey(name) ) {
                params[name] = ParamsHelper.resolveFromCode(decl, configParams[name])
            }
            else if( decl.defaultValue != null ) {
                params[name] = ParamsHelper.resolveFromCode(decl, decl.defaultValue)
            }
            else if( decl.optional ) {
                params[name] = null
            }
            else {
                throw new ScriptRuntimeException("Parameter `$name` is required but was not specified on the command line, params file, or config")
            }

            final expectedType = TypeHelper.getRawType(decl.type)
            final actualType = params[name]?.getClass()
            if( actualType != null && !ParamsHelper.isAssignableFrom(expectedType, actualType) )
                throw new ScriptRuntimeException("Parameter `$name` with type ${Types.getName(decl.type)} cannot be assigned to ${params[name]} [${Types.getName(actualType)}]")
        }

        // propagate resolved params to all scripts for legacy compatibility
        if( !session.binding.getScriptPath() )
            session.binding.setParams(params, true)
        for( final scriptPath : ScriptMeta.allScriptNames().values() ) {
            if( !scriptPath )
                continue
            final script = ScriptMeta.getScriptByPath(scriptPath)
            script.binding.setParams(params, true)
        }
    }

}
