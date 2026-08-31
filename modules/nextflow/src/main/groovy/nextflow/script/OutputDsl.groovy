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
import groovyx.gpars.dataflow.DataflowVariable
import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import nextflow.extension.CH
import nextflow.extension.DumpHelper
import nextflow.extension.PublishOp
/**
 * Implements the DSL for publishing workflow outputs
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class OutputDsl {

    private BaseScript owner

    private ScriptBinding.ParamsMap params

    private Map<String,Map> declarations = [:]

    OutputDsl(BaseScript owner=null, ScriptBinding.ParamsMap params=null) {
        this.owner = owner
        this.params = params
    }

    /**
     * The output directives of a pipeline can refer to its params, so the
     * params of the pipeline are resolved against the output block rather
     * than against the script that declares it, which is shared by every
     * execution of the pipeline.
     *
     * @param name
     */
    @Override
    Object getProperty(String name) {
        if( name == 'params' && params != null )
            return params
        return super.getProperty(name)
    }

    private Map<String,DataflowVariable> dataflowOutputs = [:]

    void declare(String name, Closure closure) {
        declare(name, null, closure)
    }

    void declare(String name, String typeName, Closure closure) {
        final options = options(closure)
        final outputType = outputType(typeName)

        // an output declared with an included output type is equivalent to
        // re-declaring each output of the included pipeline
        if( outputType ) {
            if( options )
                throw new ScriptRuntimeException("Workflow output '${name}' with type ${typeName} cannot declare its own output directives -- declare each output separately instead")
            for( final entry : outputType.getDeclarations() )
                declare0(entry.key, entry.value)
            return
        }

        declare0(name, options)
    }

    private void declare0(String name, Map options) {
        if( declarations.containsKey(name) )
            throw new ScriptRuntimeException("Workflow output '${name}' is declared more than once in the workflow output block")

        declarations[name] = options
    }

    private static Map options(Closure closure) {
        final dsl = new DeclareDsl()
        final cl = (Closure)closure.clone()
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        cl.setDelegate(dsl)
        cl.call()

        return dsl.getOptions()
    }

    private OutputTypeDef outputType(String typeName) {
        if( !typeName || !owner )
            return null
        final component = ScriptMeta.get(owner).getComponent(typeName)
        return component instanceof OutputTypeDef ? component : null
    }

    Map<String,Map> getDeclarations() { declarations }

    void apply(Session session) {
        final outputs = session.outputs
        final defaults = session.config.navigate('workflow.output', Collections.emptyMap()) as Map

        // make sure every output was assigned
        for( final name : declarations.keySet() ) {
            if( !outputs.containsKey(name) )
                throw new ScriptRuntimeException("Workflow output '${name}' was declared in the output block but not assigned in the workflow")
        }

        for( final name : outputs.keySet() ) {
            if( !declarations.containsKey(name) )
                throw new ScriptRuntimeException("Workflow output '${name}' was assigned in the workflow but not declared in the output block")
        }

        // create publish op for each output
        for( final name : outputs.keySet() ) {
            final source = outputs[name]
            final overrides = declarations[name] ?: Collections.emptyMap()
            final opts = publishOptions(name, defaults, overrides)

            if( opts.enabled == null || opts.enabled )
                dataflowOutputs[name] = new PublishOp(session, name, CH.getReadChannel(source), opts).apply()
        }

        // print workflow outputs on run completion
        session.workflowMetadata.onComplete {
            if( !session.isSuccess() )
                return
            final output = getOutput()
            if( session.outputFormat == 'json' )
                session.printConsole(DumpHelper.prettyPrintJson(output), true)
            else if( session.outputFormat == 'text' )
                printOutput(session, output)
        }
    }

    private Map publishOptions(String name, Map defaults, Map overrides) {
        final opts = defaults + overrides
        if( opts.containsKey('ignoreErrors') )
            opts.failOnError = !opts.remove('ignoreErrors')
        if( !opts.containsKey('overwrite') )
            opts.overwrite = 'standard'

        final path = opts.path as String ?: '.'
        if( path.startsWith('/') )
            throw new ScriptRuntimeException("Invalid path '${path}' for workflow output '${name}' -- it must be a relative path")
        opts.path = path

        if( opts.index && !(opts.index as Map).path )
            throw new ScriptRuntimeException("Index file definition for workflow output '${name}' is missing `path` option")

        return opts
    }

    private static void printOutput(Session session, Map<String,Object> output) {
        if( output.size() == 1 && output.keySet().first() == '$out' ) {
            session.printConsole(output.values().first().toString())
            return
        }
        final outputDir = session.outputDir?.toUriString()
        final sb = new StringBuilder()
        sb.append('\n')
        sb.append("Outputs:\n")
        if( outputDir ) {
            sb.append('\n')
            sb.append("  ${outputDir}\n")
        }
        for( final outputName : output.keySet() ) {
            final outputValue = output[outputName]
            sb.append('\n')
            if( outputValue instanceof Collection ) {
                final items = outputValue as List
                final maxItems = items.size() > 20 ? 10 : items.size()
                sb.append("  ${outputName}:\n")
                for( final item : items.subList(0, maxItems) )
                    sb.append("    - ${normalizeOutput(item, outputDir)}")
                if( maxItems < items.size() )
                    sb.append("    - ... (${items.size() - maxItems} more items)\n")
            }
            else if( outputValue instanceof Map ) {
                sb.append("  ${outputName}:\n")
                for( final mapEntry : outputValue.entrySet() )
                    sb.append("    ${mapEntry.key}: ${normalizeOutput(mapEntry.value, outputDir)}")
            }
            else {
                sb.append("  ${outputName}: ${normalizeOutput(outputValue, outputDir)}")
            }
        }
        session.printConsole(sb.toString())
    }

    private static String normalizeOutput(Object value, String outputDir) {
        final result = DumpHelper.prettyPrintYaml(value, style: 'flow')
        return outputDir
            ? result.replace(outputDir + '/', '')
            : result
    }

    Map<String,Object> getOutput() {
        dataflowOutputs.collectEntries { name, dv -> [name, dv.get()] }
    }

    static class DeclareDsl {

        private Map opts = [:]

        void contentType(String value) {
            setOption('contentType', value)
        }

        void contentType(boolean value) {
            setOption('contentType', value)
        }

        void enabled(boolean value) {
            setOption('enabled', value)
        }

        void ignoreErrors(boolean value) {
            setOption('ignoreErrors', value)
        }

        void index(Closure closure) {
            final dsl = new IndexDsl()
            final cl = (Closure)closure.clone()
            cl.setResolveStrategy(Closure.DELEGATE_FIRST)
            cl.setDelegate(dsl)
            cl.call()
            setOption('index', dsl.getOptions())
        }

        void label(CharSequence value) {
            final opts = getOptions()
            final current = opts.get('labels')
            if( current instanceof List )
                current.add(value)
            else
                opts.put('labels', [value])
        }

        void mode(String value) {
            setOption('mode', value)
        }

        void overwrite(boolean value) {
            setOption('overwrite', value)
        }

        void overwrite(String value) {
            setOption('overwrite', value)
        }

        void path(String value) {
            setOption('path', value)
        }

        void path(Closure value) {
            setOption('path', '.')
            setOption('pathResolver', value)
        }

        void storageClass(String value) {
            setOption('storageClass', value)
        }

        void tags(Map value) {
            setOption('tags', value)
        }

        private void setOption(String name, Object value) {
            if( opts.containsKey(name) )
                throw new ScriptRuntimeException("Publish option `${name}` cannot be defined more than once for a workflow output")
            opts[name] = value
        }

        Map getOptions() {
            return opts
        }

    }

    static class IndexDsl {

        private Map opts = [:]

        void header(boolean value) {
            setOption('header', value)
        }

        void header(List<String> value) {
            setOption('header', value)
        }

        void header(String... value) {
            setOption('header', value as List)
        }

        void path(String value) {
            setOption('path', value)
        }

        void sep(String value) {
            setOption('sep', value)
        }

        private void setOption(String name, Object value) {
            if( opts.containsKey(name) )
                throw new ScriptRuntimeException("Index option `${name}` cannot be defined more than once for a given index definition")
            opts[name] = value
        }

        Map getOptions() {
            return opts
        }

    }

}
