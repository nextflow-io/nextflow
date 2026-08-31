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
import nextflow.extension.Bolts
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

    private BaseScript owner

    private Class clazz

    private Map<String,Param> declarations = [:]

    ParamsDsl(BaseScript owner, Class clazz) {
        this.owner = owner
        this.clazz = clazz
    }

    void declare(String name, boolean optional, Object defaultValue = null) {
        final type = clazz.getField(name).getGenericType()
        if( defaultValue == null && type == Boolean )
            defaultValue = false
        declarations[name] = new Param(name, type, optional, defaultValue)
    }

    Map<String,Param> getDeclarations() { declarations }

    /**
     * Resolve the declared params against the command line and config,
     * and return them.
     *
     * @param session
     */
    ScriptBinding.ParamsMap apply(Session session) {
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
                // a nested param (e.g. `--rnaseq.aligner`) overrides only the
                // fields it names, keeping the rest of the config value
                final value = cliParams[name] instanceof Map && configParams[name] instanceof Map
                    ? Bolts.deepMerge((Map)configParams[name], (Map)cliParams[name])
                    : cliParams[name]
                params[name] = ParamsHelper.resolveFromCli(decl, value)
            }
            else if( configParams.containsKey(name) ) {
                params[name] = ParamsHelper.resolveFromCode(decl, configParams[name])
            }
            else if( decl.defaultValue != null ) {
                params[name] = ParamsHelper.resolveFromCode(decl, decl.defaultValue)
            }
            else {
                params[name] = ParamsHelper.emptyRecord(decl)
            }

            if( params[name] == null && !decl.optional ) {
                throw new ScriptRuntimeException("Parameter `$name` is required but no value was provided")
            }

            ParamsHelper.checkAssignable(decl, params[name])
        }

        // propagate resolved params to all scripts for legacy compatibility
        if( !session.binding.getScriptPath() )
            session.binding.setParams(params, true)
        for( final scriptPath : ScriptMeta.allScriptNames().values() ) {
            if( !scriptPath )
                continue
            final script = ScriptMeta.getScriptByPath(scriptPath)
            // a script that declares its own params block resolves its params
            // independently -- a pipeline included by this script receives its
            // params from the pipeline call, not from the calling pipeline
            if( !script.is(owner) && script.getParamDeclarations() )
                continue
            script.binding.setParams(params, true)
        }

        // the resolved params of the owner script also include any config
        // param that the script does not declare
        return owner.binding.getParams()
    }

}
