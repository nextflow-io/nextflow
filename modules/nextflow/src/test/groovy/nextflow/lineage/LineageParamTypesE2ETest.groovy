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
package nextflow.lineage

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Session
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.model.v1beta1.TaskOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptFile
import nextflow.script.ScriptLoaderFactory
import spock.lang.TempDir
import spock.lang.Timeout
import test.Dsl2Spec

/**
 * End-to-end test of the parameter types recorded in lineage: typed processes
 * record their declared types, legacy processes record their qualifiers.
 */
@Timeout(120)
class LineageParamTypesE2ETest extends Dsl2Spec {

    @TempDir
    Path tempDir

    def 'should record declared types for typed processes'() {
        given:
        final ref = tempDir.resolve('ref.txt'); ref.text = 'ref'
        final script = tempDir.resolve('main.nf')
        script.text = """
            nextflow.enable.types = true

            record Sample {
                id: String
                reads: Path
            }

            process MAKE {
                input:
                id: String
                n: Integer
                meta: Map<String,String>
                ref: Path
                opt: Path?

                output:
                sample: Sample = record(id: id, reads: file('a.txt'))
                pair: Tuple<String,Path> = tuple(id, file('a.txt'))
                all: Set<Path> = files('*.txt')
                missing: Path? = file('missing.txt', optional: true)
                out: String = stdout()

                script:
                "touch a.txt b.txt; echo hi"
            }

            process USE_RECORD {
                input:
                sample: Sample

                output:
                file('c.txt')

                script:
                "touch c.txt"
            }

            process USE_TUPLE {
                input:
                tuple(id: String, f: Path)

                output:
                file('d.txt')

                script:
                "touch d.txt"
            }

            process USE_LIST {
                input:
                files: List<Path>

                output:
                files

                script:
                "true"
            }

            workflow {
                make = MAKE(channel.of('a'), 1, [k: 'v'], file('${ref}'), null)
                USE_RECORD(make.sample)
                USE_TUPLE(make.pair)
                USE_LIST(make.all.map { fs -> fs.toList() })
            }
            """

        when:
        final store = runWithLineage(script)

        then:
        inputTypes(store, 'MAKE') == [id: 'String', n: 'Integer', meta: 'Map', ref: 'Path', opt: 'Path']
        outputTypes(store, 'MAKE') == [sample: 'Record', pair: 'Tuple', all: 'Set', missing: 'Path', out: 'String']
        and:
        inputTypes(store, 'USE_RECORD') == [sample: 'Record']
        inputTypes(store, 'USE_TUPLE') == [id: 'String', f: 'Path']
        inputTypes(store, 'USE_LIST') == [files: 'List']
        and: 'an unnamed output is not declared, so it records the value type'
        outputTypes(store, 'USE_RECORD') == ['$out': 'Path']
        and: 'a re-emitted input keeps its declared type'
        outputTypes(store, 'USE_LIST') == [files: 'List']
    }

    def 'should record qualifiers for legacy processes'() {
        given:
        final ref = tempDir.resolve('ref.txt'); ref.text = 'ref'
        final other = tempDir.resolve('other.txt'); other.text = 'other'
        final script = tempDir.resolve('main.nf')
        script.text = """
            process LEGACY {
                input:
                val id
                path ref
                env 'FOO'
                stdin
                tuple val(key), path(f)

                output:
                val id
                path 'a.txt'
                env 'BAR'
                stdout
                eval 'echo hi'
                tuple val(key), path('a.txt')

                script:
                "cat - > /dev/null; touch a.txt; export BAR=baz; echo out"
            }

            workflow {
                LEGACY('a', file('${ref}'), 'foo', 'in', channel.of(['x', file('${other}')]))
            }
            """

        when:
        final store = runWithLineage(script)

        then:
        loadTaskRun(store, 'LEGACY').input.collect { p -> [p.name, p.type] } == [
            ['id', 'val'], ['ref', 'path'], ['FOO', 'env'], ['-', 'stdin'], ['key', 'val'], ['f', 'path']
        ]
        loadTaskOutput(store, 'LEGACY').output*.type == ['val', 'path', 'env', 'stdout', 'eval', 'val', 'path']
    }

    private DefaultLinStore runWithLineage(Path script) {
        final config = [
            workDir: tempDir.resolve('work').toString(),
            lineage: [enabled: true, store: [location: tempDir.resolve('lineage').toString()]] ]

        final session = new Session(config)
        session.setBinding(new ScriptBinding())
        session.init(new ScriptFile(script), null, null, null)

        final store = new DefaultLinStore()
        store.open(LineageConfig.create(session))
        final observers = Session.getDeclaredField('observersV2')
        observers.setAccessible(true)
        observers.set(session, (observers.get(session) as List) + [new LinObserver(session, store)])

        session.start()
        final loader = ScriptLoaderFactory.create(session)
        loader.parse(script)
        loader.runScript()
        session.fireDataflowNetwork()
        session.await()
        session.destroy()
        if( session.error )
            throw session.error
        return store
    }

    private String taskKey(DefaultLinStore store, String process) {
        final keys = Files.list(tempDir.resolve('lineage')).toList()
            .collect { it.fileName.toString() }
            .findAll { key ->
                final record = store.load(key)
                record instanceof TaskRun && record.name.startsWith(process + ' ')
            }
        assert keys.size() == 1, "Expected exactly one task for process ${process}, found ${keys.size()}"
        return keys[0]
    }

    private TaskRun loadTaskRun(DefaultLinStore store, String process) {
        return store.load(taskKey(store, process)) as TaskRun
    }

    private TaskOutput loadTaskOutput(DefaultLinStore store, String process) {
        return store.load(taskKey(store, process) + '#output') as TaskOutput
    }

    private Map<String,String> inputTypes(DefaultLinStore store, String process) {
        return loadTaskRun(store, process).input.collectEntries { p -> [p.name, p.type] }
    }

    private Map<String,String> outputTypes(DefaultLinStore store, String process) {
        return loadTaskOutput(store, process).output.collectEntries { p -> [p.name, p.type] }
    }
}
