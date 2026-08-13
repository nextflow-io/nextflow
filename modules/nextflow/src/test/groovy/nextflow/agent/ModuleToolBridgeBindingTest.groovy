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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Nextflow
import nextflow.Session
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpecFactory
import nextflow.script.ProcessDef
import nextflow.script.ProcessEntryHandler
import nextflow.script.ScriptMeta
import nextflow.script.parser.v2.ScriptLoaderV2
import test.Dsl2Spec

/**
 * Asserts the agent tool bridge marshals the LLM's tool-call args into the module's input
 * channel using the SAME logic as {@code nextflow module run} -- i.e. through
 * {@link ProcessEntryHandler#getProcessArguments(ProcessDef, Map, ModuleSpec)} -- so the
 * {@code meta} map, file coercion and tuple assembly are identical.
 *
 * For a classic-DSL2 {@code tuple val(meta), path(reads)} input (described by a sibling
 * {@code meta.yml}), the flattened args {@code {meta:[id:'s1'], reads:'/abs/x'}} must build the
 * channel value {@code [[id:'s1'], file('/abs/x')]} -- exactly what the bridge binds onto the
 * pre-wired input queue.
 */
class ModuleToolBridgeBindingTest extends Dsl2Spec {

    /** Compile a classic-DSL2 (V1) module to its single ProcessDef, like AgentDef does. */
    private ProcessDef loadModuleProcess(Path modPath) {
        final session = new Session()
        final loader = new ScriptLoaderV2(session)
        loader.setModule(true)
        loader.parse(modPath)
        loader.runScript()
        final meta = ScriptMeta.get(loader.getScript())
        return meta.getProcess(meta.getProcessNames().first())
    }

    def 'should build the same channel value as module run for a tuple val(meta), path(reads) input'() {
        given:
        final dir = Files.createTempDirectory('test')
        final modPath = dir.resolve('main.nf')
        modPath.text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("out.txt"), emit: report

                script:
                """
                cat ${reads} > out.txt
                """
            }
            '''.stripIndent()

        and: 'a sibling meta.yml describing the tuple I/O (map + file)'
        final metaPath = dir.resolve('meta.yml')
        metaPath.text = '''\
            name: echo_tool
            description: Copy the input reads to an output file
            input:
              - - name: meta
                  type: map
                  description: sample meta
                - name: reads
                  type: file
                  description: input reads
            output:
              - - name: meta
                  type: map
                  description: sample meta
                - name: outfile
                  type: file
                  description: the output
            '''.stripIndent()

        and:
        final ProcessDef proc = loadModuleProcess(modPath)
        final ModuleSpec spec = ModuleSpecFactory.fromYaml(metaPath)

        when: 'the bridge marshals the flattened LLM args via the module-run binding'
        final args = ProcessEntryHandler.getProcessArguments(proc, [meta: [id: 's1'], reads: '/abs/x'], spec)

        then: 'one element per input channel: the tuple is assembled to [meta-map, file(reads)]'
        args.size() == 1
        args[0] instanceof List
        (args[0] as List).size() == 2
        (args[0] as List)[0] == [id: 's1']
        (args[0] as List)[1] == Nextflow.file('/abs/x')
        (args[0] as List)[1] instanceof Path
    }

    def 'should raise IllegalArgumentException for a missing required arg (LLM-recoverable contract)'() {
        given:
        final dir = Files.createTempDirectory('test')
        final modPath = dir.resolve('main.nf')
        modPath.text = '''
            process echo_tool {
                input:
                tuple val(meta), path(reads)

                output:
                tuple val(meta), path("out.txt"), emit: report

                script:
                """
                cat ${reads} > out.txt
                """
            }
            '''.stripIndent()

        and:
        final metaPath = dir.resolve('meta.yml')
        metaPath.text = '''\
            name: echo_tool
            input:
              - - name: meta
                  type: map
                - name: reads
                  type: file
            '''.stripIndent()

        and:
        final ProcessDef proc = loadModuleProcess(modPath)
        final ModuleSpec spec = ModuleSpecFactory.fromYaml(metaPath)

        when: 'a required `meta` arg is missing'
        ProcessEntryHandler.getProcessArguments(proc, [reads: '/abs/x'], spec)

        then: 'the binding throws -- the bridge dispatcher turns this into a {"error":...} tool result'
        thrown(IllegalArgumentException)
    }
}
