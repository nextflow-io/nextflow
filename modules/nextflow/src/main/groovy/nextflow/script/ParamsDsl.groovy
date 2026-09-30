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

    Map<String,Param> getDeclarations() { declarations }

    void apply(Session session) {
        final params = ParamsHelper.resolveParams(declarations.values(), session.cliParams ?: [:], session.configParams ?: [:])

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
