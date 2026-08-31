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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
/**
 * Models the workflow output definition
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class OutputDef {

    private BaseScript owner

    private Closure closure

    OutputDef(BaseScript owner, Closure closure) {
        this.owner = owner
        this.closure = closure
    }

    void apply(Session session, ScriptBinding.ParamsMap params) {
        dsl(params).apply(session)
    }

    /**
     * The declared outputs, keyed by name, with their output directives
     * resolved against the given params.
     *
     * The output block is evaluated on each use, because the output
     * directives can refer to the params of the pipeline, and a pipeline
     * can be executed more than once in a run.
     *
     * @param params
     */
    Map<String,Map> getDeclarations(ScriptBinding.ParamsMap params) {
        return dsl(params).getDeclarations()
    }

    private OutputDsl dsl(ScriptBinding.ParamsMap params) {
        final dsl = new OutputDsl(owner, params)
        final cl = (Closure)closure.clone()
        cl.setDelegate(dsl)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.call()
        return dsl
    }

}
