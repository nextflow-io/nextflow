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

    private Map<String,Optional<?>> declarations = [:]

    void declare(String name) {
        declare0(name, Optional.empty())
    }

    void declare(String name, Object defaultValue) {
        declare0(name, Optional.of(defaultValue))
    }

    private void declare0(String name, Optional<?> defaultValue) {
        if( declarations.containsKey(name) )
            throw new ScriptRuntimeException("Parameter '${name}' is declared more than once in the workflow params definition")

        declarations[name] = defaultValue
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
            final defaultValue = declarations[name]
            if( cliParams.containsKey(name) )
                params[name] = cliParams[name]
            else if( configParams.containsKey(name) )
                params[name] = configParams[name]
            else if( defaultValue.isPresent() )
                params[name] = defaultValue.get()
            else
                throw new ScriptRuntimeException("Parameter `$name` is required but was not specified on the command line, params file, or config")
        }

        session.binding.setParams(params)
    }

}
