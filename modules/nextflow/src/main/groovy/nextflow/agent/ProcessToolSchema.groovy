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
package nextflow.agent

import groovy.transform.CompileStatic
import nextflow.script.ProcessConfigV2
import nextflow.script.ProcessDef
import nextflow.script.params.v2.ProcessInput
import nextflow.script.params.v2.ProcessOutput
import nextflow.script.params.v2.ProcessTupleInput

/**
 * Derives portable JSON-schema {@link Map}s describing the declared typed inputs
 * and outputs of a {@link ProcessDef}, so an in-scope process can be exposed to
 * the LLM as a tool. The resulting maps share the shape produced by
 * {@link RecordSchema#of} (e.g.
 * {@code [type:'object', properties:[name:[type:'string']], required:['name'], additionalProperties:false]}).
 *
 * For Phase 2 only SCALAR inputs/outputs are supported (String, integer, number,
 * boolean); any other declared kind (tuple, path/file, map/meta, record) raises
 * a loud {@link IllegalArgumentException} so unsupported processes fail fast
 * rather than silently producing a wrong schema.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ProcessToolSchema {

    /**
     * Build the JSON-schema map for the process inputs. Each declared input
     * contributes a {@code name -> {type: ...}} property; non-optional inputs are
     * listed under {@code required}.
     */
    static Map inputSchema(ProcessDef proc) {
        final config = configOf(proc)
        final properties = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()

        for( final param : config.getInputs().getParams() ) {
            if( param instanceof ProcessTupleInput )
                throw unsupported(proc.getName(), 'input', '(tuple)', 'Tuple')

            final name = param.getName()
            final type = param.getType()
            final fragment = RecordSchema.scalarFragment(type)
            if( fragment == null )
                throw unsupported(proc.getName(), 'input', name, type)

            properties.put(name, fragment)
            if( !((ProcessInput) param).isOptional() )
                required.add(name)
        }

        return assemble(properties, required)
    }

    /**
     * Build the JSON-schema map for the process outputs. For Phase 2 a single
     * scalar output is supported; a bare typed value with no declared name is
     * exposed under the key {@code result}.
     */
    static Map outputSchema(ProcessDef proc) {
        final config = configOf(proc)
        final properties = new LinkedHashMap<String,Object>()
        final required = new ArrayList<String>()

        for( final param : config.getOutputs().getParams() ) {
            final name = outputKey(param)
            final type = param.getType()
            final fragment = RecordSchema.scalarFragment(type)
            if( fragment == null )
                throw unsupported(proc.getName(), 'output', name, type)

            properties.put(name, fragment)
            required.add(name)
        }

        return assemble(properties, required)
    }

    private static ProcessConfigV2 configOf(ProcessDef proc) {
        final config = proc.getProcessConfig()
        if( !(config instanceof ProcessConfigV2) )
            throw new IllegalArgumentException("Tool `${proc.getName()}` must be a typed process to be used as an agent tool")
        return (ProcessConfigV2) config
    }

    /**
     * The key used for an output property. A named output (e.g. {@code result: String})
     * uses its declared name; a bare typed value lowers to the synthetic name
     * {@code $out}, for which the key {@code result} is used instead.
     */
    private static String outputKey(ProcessOutput param) {
        final name = param.getName()
        return ( name == null || name == '$out' ) ? 'result' : name
    }

    private static Map assemble(Map<String,Object> properties, List<String> required) {
        final result = new LinkedHashMap<String,Object>()
        result.put('type', 'object')
        result.put('properties', properties)
        result.put('required', required)
        result.put('additionalProperties', false)
        return result
    }

    private static IllegalArgumentException unsupported(String proc, String kind, String name, Object type) {
        final typeName = type instanceof Class ? ((Class) type).getSimpleName() : type
        return new IllegalArgumentException("Tool `${proc}` ${kind} `${name}` of type ${typeName} is not yet supported as an agent tool (Phase 3)")
    }

}
