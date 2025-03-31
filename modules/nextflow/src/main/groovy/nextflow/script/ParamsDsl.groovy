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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
/**
 * Implements the DSL for defining workflow params
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ParamsDsl {

    private Map<String,Map> declarations = [:]

    void declare(String name, Closure closure) {
        if( declarations.containsKey(name) )
            throw new ScriptRuntimeException("Parameter '${name}' is declared more than once in the workflow params definition")

        final dsl = new DeclareDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()

        declarations[name] = dsl.getOptions()
    }

    void apply(Session session) {
        final cliParams = session.params
        final configParams = session.config.params as Map

        for( final name : cliParams.keySet() ) {
            // TODO: doesn't work because CLI params are added to config params
            // TODO: need to use config params before CLI params were added
            if( !declarations.containsKey(name) && !configParams.containsKey(name) )
                throw new ScriptRuntimeException("Parameter `$name` was specified on the command line or params file but is not declared in the script or config")
        }

        final params = new HashMap<String,?>()
        for( final name : declarations.keySet() ) {
            final decl = declarations[name]
            if( configParams.containsKey(name) )
                params[name] = configParams[name]
            else if( decl.containsKey('defaultValue') )
                params[name] = decl.defaultValue
            else
                throw new ScriptRuntimeException("Parameter `$name` is required but was not specified on the command line, params file, or config")
        }

        session.binding.setParams(params)
    }

    static class DeclareDsl {

        private Map opts = [:]

        void defaultValue(Object value) {
            setOption('defaultValue', value)
        }

        void description(String value) {
            setOption('description', value)
        }

        private void setOption(String name, Object value) {
            if( opts.containsKey(name) )
                throw new ScriptRuntimeException("Parameter directive `${name}` cannot be defined more than once for a given parameter")
            opts[name] = value
        }

        Map getOptions() {
            return opts
        }

    }

}
