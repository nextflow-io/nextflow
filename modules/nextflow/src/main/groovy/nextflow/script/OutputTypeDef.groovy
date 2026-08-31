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
import nextflow.exception.ScriptRuntimeException

/**
 * Models the output block of a pipeline included into another script
 * as a record type:
 *
 *   include { output as RnaseqOutput } from './pipelines/rnaseq.nf'
 *
 * Declaring an output with this type is equivalent to re-declaring each
 * output of the included pipeline, along with its output directives.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class OutputTypeDef extends ComponentDef {

    private BaseScript pipeline

    private String name

    private Map<String,Map> declarations

    OutputTypeDef(BaseScript pipeline, String name='output') {
        this.pipeline = pipeline
        this.name = name
    }

    BaseScript getPipeline() { pipeline }

    String getName() { name }

    String getType() { 'output type' }

    OutputTypeDef cloneWithName(String name) {
        final result = (OutputTypeDef)clone()
        result.@name = name
        return result
    }

    /**
     * The outputs of the included pipeline, keyed by name, with their
     * output directives resolved against the params of the pipeline call.
     *
     * The declarations are supplied by the pipeline when it is executed,
     * because the output directives can refer to its params.
     *
     * @param declarations
     */
    void setDeclarations(Map<String,Map> declarations) {
        this.declarations = declarations
    }

    Map<String,Map> getDeclarations() {
        if( declarations == null )
            throw new ScriptRuntimeException("Output type `${name}` cannot be resolved because the corresponding pipeline was not invoked")
        return declarations
    }

    @Override
    Object invoke_a(Object[] args) {
        throw new ScriptRuntimeException("Output type `${name}` cannot be invoked -- it can only be used to declare a workflow output")
    }

}
