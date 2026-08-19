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

import nextflow.Session
import nextflow.executor.SimpleFileCopyStrategy
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ScriptBinding
import nextflow.script.ScriptLoaderFactory
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import spock.lang.Specification
import spock.lang.TempDir
import spock.lang.Timeout

/**
 * Parity of an agent's typed `path` I/O with a typed process's: a declared `Path` input is staged
 * into the task work dir (and materialized there for an in-JVM agent), and a `path` output is
 * collected out of it.
 *
 * Every test drives a REAL {@link Session} — real {@code ExecutorFactory}, real
 * {@link nextflow.executor.local.AgentExecutor}, real work dirs — because the mock harness never
 * instantiates the agent executor and never stages anything. The LLM is a stub runner, so no test
 * here makes a model call.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(60)
class AgentPathIoTest extends Specification {

    @TempDir
    Path tempDir

    def setup() {
        // `errorShown` is a static one-shot: without this, the FIRST test in the class that
        // aborts consumes it and every later failure returns a bare ErrorStrategy instead of a
        // TaskFault, so the session never aborts and a negative test cannot observe the failure
        TaskProcessor.reset()
    }

    def cleanup() {
        AgentRunnerProvider.testRunner = null
    }

    // -----------------------------------------------------------------------
    // inputs
    // -----------------------------------------------------------------------

    def 'should stage a typed Path input into the agent task work dir'() {
        given:
        final input = tempDir.resolve('contigs.fa')
        input.text = '>chr1\nACGT\n'
        and:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, """
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path
                output:
                answer: String
                prompt: "Inspect \${contigs}"
            }

            workflow {
                qa(channel.of(file('${input}')))
            }
            """)

        then: 'the declaration registered a file input, exactly as a typed process does'
        final task = tasks.completed[0]
        task.inputFiles.size() == 1
        task.inputFiles[0].stageName == 'contigs.fa'
        Files.isSameFile(task.inputFiles[0].sourcePath, input)
        task.getInputFilesMap().keySet() == ['contigs.fa'] as Set

        and: 'an in-JVM agent has no wrapper script, so the handler materialized it into the work dir'
        final staged = task.workDir.resolve('contigs.fa')
        Files.exists(staged)
        Files.isSymbolicLink(staged)
        Files.isSameFile(staged, input)

        and: 'the model was given the work-dir-relative name, NOT a driver-side absolute path'
        captured.inputJson == '"contigs.fa"'
        and: 'and the prompt interpolation agrees with it -- one input, one rendering'
        captured.prompt == 'Inspect contigs.fa'
    }

    def 'should let an in-JVM agent OPEN its staged input through the fs: tools'() {
        given: 'the staging is a symlink, and SandboxGuard resolves symlinks before testing containment'
        // a text-like suffix so `read` inlines the content rather than returning an opaque handle
        final input = tempDir.resolve('contigs.txt')
        input.text = '>chr1\nACGT\n'
        and:
        String readResult = null
        String listResult = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            readResult = req.dispatch.call('read', '{"path":"contigs.txt"}')
            listResult = req.dispatch.call('ls', '{"path":"."}')
            'ok'
        } as AgentRunner

        when:
        runWithObserver(taskProbe().probe, """
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                tools 'fs:read', 'fs:ls'
                input:
                contigs: Path
                output:
                answer: String
                prompt: "Inspect \${contigs}"
            }

            workflow {
                qa(channel.of(file('${input}')))
            }
            """)

        then: 'the name the agent was handed is a name it can actually open'
        !readResult.contains('outside sandbox')
        readResult.contains('ACGT')
        and: 'and one it can SEE -- not an opaque `link` entry with no size'
        listResult.contains('"name":"contigs.txt"')
        listResult.contains('"type":"file"')
    }

    def 'should stage an agent Path input under the same name a process stages it'() {
        given:
        final input = tempDir.resolve('reads_1.fq')
        input.text = 'x'
        and:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()

        when: 'an agent and a process declare the very same input'
        runWithObserver(tasks.probe, """
            nextflow.enable.types = true

            process CHECK {
                input:
                contigs: Path
                output:
                answer: String = 'x'
                script:
                'true'
            }

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                CHECK(channel.of(file('${input}')))
                qa(channel.of(file('${input}')))
            }
            """)

        then: 'both stage it, and both stage it under the same name'
        final names = tasks.completed.collectEntries { [(it.processor.name): it.getStagedInputs()] }
        names['CHECK'] == ['reads_1.fq']
        names['qa'] == ['reads_1.fq']
    }

    def 'should produce the same staged-input shape as a process with the same declaration'() {
        given: 'one input of every shape the stager inference recognises, plus one it must ignore'
        final contigs = tempDir.resolve('contigs.fa'); contigs.text = 'c'
        final r1 = tempDir.resolve('r1.fq'); r1.text = '1'
        final r2 = tempDir.resolve('r2.fq'); r2.text = '2'
        final e1 = tempDir.resolve('e1.txt'); e1.text = 'a'
        final e2 = tempDir.resolve('e2.txt'); e2.text = 'b'
        and:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()
        // FileHolder carries the whole staging decision -- source object, store path and stage
        // name -- and is @EqualsAndHashCode, so comparing the sorted holders compares the SHAPE,
        // not merely the names. Sorted because declaration order is not part of the contract.
        final shape = { TaskRun t -> t.inputFiles.toSorted { it.stageName } }

        when: 'an agent and a process declare the very same typed inputs'
        runWithObserver(tasks.probe, """
            nextflow.enable.types = true

            record Pair {
                id: String
                r1: Path
                r2: Path
            }

            process CHECK {
                input:
                contigs: Path
                pair: Pair
                extras: List<Path>
                tag: String
                output:
                answer: String = 'x'
                script:
                'true'
            }

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path
                pair: Pair
                extras: List<Path>
                tag: String
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                CHECK(
                    channel.value(file('${contigs}')),
                    channel.value(record(id: 's1', r1: file('${r1}'), r2: file('${r2}'))),
                    channel.value([file('${e1}'), file('${e2}')]),
                    channel.value('t') )
                qa(
                    channel.value(file('${contigs}')),
                    channel.value(record(id: 's1', r1: file('${r1}'), r2: file('${r2}'))),
                    channel.value([file('${e1}'), file('${e2}')]),
                    channel.value('t') )
            }
            """)

        then: 'the agent stages exactly what the process stages -- same files, same stage names'
        final byName = tasks.completed.collectEntries { [(it.processor.name): it] }
        shape(byName['qa']) == shape(byName['CHECK'])

        and: 'and the equality is not vacuous: the scalar, both record fields and both collection'
        // elements are staged, while the non-Path input contributes nothing
        shape(byName['CHECK'])*.stageName == ['contigs.fa', 'e1.txt', 'e2.txt', 'r1.fq', 'r2.fq']
    }

    def 'should stage every Path field of a record input'() {
        given:
        final fa = tempDir.resolve('sample.fa')
        fa.text = 'a'
        final idx = tempDir.resolve('sample.idx')
        idx.text = 'b'
        and:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, """
            nextflow.enable.types = true

            record Sample {
                id: String
                seq: Path
                index: Path
            }

            agent qa {
                model 'openai/gpt-4o'
                input:
                sample: Sample
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                qa(channel.of(record(id: 's1', seq: file('${fa}'), index: file('${idx}'))))
            }
            """)

        then: 'the record recursion registered one stager per Path field'
        final task = tasks.completed[0]
        task.getStagedInputs().toSet() == ['sample.fa', 'sample.idx'] as Set

        and: 'the model sees the staged names inside the record, not driver-side paths'
        captured.inputJson == '{"id":"s1","seq":"sample.fa","index":"sample.idx"}'
    }

    def 'should admit a null value for an optional Path input and stage nothing'() {
        given:
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> captured = req; 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path?
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                qa(channel.of(null))
            }
            ''')

        then: 'the task is constructed, and the absent input stages nothing'
        tasks.completed.size() == 1
        tasks.completed[0].inputFiles.isEmpty()
        captured.inputJson == 'null'
    }

    def 'should reject a null value for a non-optional Path input'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner

        when:
        runWithObserver(taskProbe().probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                qa(channel.of(null))
            }
            ''')

        then:
        def e = thrown(Exception)
        allMessages(e).contains('cannot be null')
    }

    // -----------------------------------------------------------------------
    // outputs
    // -----------------------------------------------------------------------

    def 'should collect a path output from the agent work dir'() {
        given: 'a runner that writes the file the prompt asked for'
        AgentRunnerRequest captured = null
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            captured = req
            Path.of(req.workDir).resolve('report.md').text = '# done\n'
            'the model chatter that must be discarded'
        } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        final session = runWithObserver(tasks.probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                report: Path = file('report.md')
                prompt: "write the report to report.md"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then: 'the unstager was registered, so the file is a real task output'
        final task = tasks.completed[0]
        task.outputFiles.size() == 1
        task.outputFiles[0] == task.workDir.resolve('report.md')
        task.outputFiles[0].text == '# done\n'

        and: 'an agent whose ONLY output is a collected file is legal: the model is given no'
        // output contract at all, and its final text is explicitly discarded
        captured.outputSchema == null
        session.error == null
    }

    def 'should surface an arity error when the agent writes a different file name'() {
        given: 'a runner that writes the WRONG name'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            Path.of(req.workDir).resolve('summary.md').text = 'oops'
            'done'
        } as AgentRunner

        when:
        runWithObserver(taskProbe().probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                report: Path = file('report.md')
                prompt: "write the report to report.md"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then:
        def e = thrown(Exception)
        // the work-dir collector rejects the missing file before the arity check is reached;
        // either way a wrong filename fails the task loudly rather than binding nothing
        allMessages(e).contains('Missing output file(s) `report.md`')
    }

    def 'should collect a files() glob output from the agent work dir'() {
        given:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            Path.of(req.workDir).resolve('a.txt').text = 'a'
            Path.of(req.workDir).resolve('b.txt').text = 'b'
            'done'
        } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                notes: Set<Path> = files('*.txt')
                prompt: "write notes as .txt files"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then:
        final task = tasks.completed[0]
        task.outputFiles.collect { it.name }.toSet() == ['a.txt', 'b.txt'] as Set
    }

    def 'should admit a missing optional path output'() {
        given: 'a runner that writes nothing at all'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'nothing to report' } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                report: Path = file(optional: true, 'report.md')
                prompt: "write report.md only if there is something to say"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''')

        then: 'the named option reached the collector, so the absent file is not a failure'
        tasks.completed.size() == 1
        tasks.completed[0].outputFiles.isEmpty()
    }

    def 'should publish a collected path output'() {
        given:
        final publishDir = Files.createTempDirectory(tempDir, 'published')
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req ->
            Path.of(req.workDir).resolve('report.md').text = '# done\n'
            'done'
        } as AgentRunner

        when:
        // D4's knock-on: a collected output is what TaskProcessor.getPublishFiles reads for a V2
        // config, so publishing an agent's file is a consequence of registering the unstager.
        // `publishDir` is not an agent DIRECTIVE (AgentBuilder.DIRECTIVES), so it can only be set
        // from the `agent` config scope -- which is the whole task ladder, minus the agent options.
        runWithObserver(taskProbe().probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                report: Path = file('report.md')
                prompt: "write the report to report.md"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''', [agent: [publishDir: [path: publishDir.toString(), mode: 'copy']]])

        then:
        publishDir.resolve('report.md').text == '# done\n'
    }

    def 'should refuse storeDir together with a collected path output on an in-JVM runner'() {
        given: 'nothing copies the work-dir file to the store dir without a wrapper script'
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'done' } as AgentRunner

        when:
        runWithObserver(taskProbe().probe, '''
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                q: String
                output:
                report: Path = file('report.md')
                prompt: "write the report to report.md"
            }

            workflow {
                qa(channel.of('hello'))
            }
            ''', [agent: [storeDir: Files.createTempDirectory(tempDir, 'store').toString()]])

        then: 'refused up front, rather than failing every run with a missing output'
        def e = thrown(Exception)
        allMessages(e).contains('cannot combine `storeDir` with a collected `path` output')
    }

    // -----------------------------------------------------------------------
    // the canonical (wrapper-script) half
    // -----------------------------------------------------------------------

    def 'should hand the staged inputs to the stage-in script a canonical agent runs'() {
        given:
        final input = tempDir.resolve('contigs.fa')
        input.text = 'c'
        and:
        AgentRunnerProvider.testRunner = { AgentRunnerRequest req -> 'ok' } as AgentRunner
        and:
        final tasks = taskProbe()

        when:
        runWithObserver(tasks.probe, """
            nextflow.enable.types = true

            agent qa {
                model 'openai/gpt-4o'
                input:
                contigs: Path
                output:
                answer: String
                prompt: "go"
            }

            workflow {
                qa(channel.of(file('${input}')))
            }
            """)

        then: 'the very inputs BashWrapperBuilder stages in and bind-mounts are populated'
        // a canonical agent differs from this one only in its BodyDef type: `task.inputFiles` is
        // filled by the task-type-agnostic V2 input resolver and consumed from the TaskBean, so
        // pinning the bean and the script it generates pins the containerized path too
        final task = tasks.completed[0]
        final bean = new TaskBean(task)
        bean.inputFiles.keySet() == ['contigs.fa'] as Set
        Files.isSameFile(bean.inputFiles['contigs.fa'], input)
        and:
        final script = new SimpleFileCopyStrategy(bean).getStageInputFilesScript(bean.inputFiles)
        script.contains("ln -sfn ${bean.inputFiles['contigs.fa']} contigs.fa".toString())
    }

    // -----------------------------------------------------------------------
    // helpers
    // -----------------------------------------------------------------------

    /** A probe capturing the completed {@link TaskRun}s, plus the list it fills. */
    private static Map taskProbe() {
        final completed = new ArrayList<TaskRun>()
        final probe = new TraceObserverV2() {
            @Override
            void onTaskComplete(TaskEvent event) { completed.add(event.handler.task) }
        }
        return [probe: probe, completed: completed]
    }

    /** Concatenate the message of a throwable and its cause chain. */
    private static String allMessages(Throwable t) {
        final sb = new StringBuilder()
        while( t != null ) {
            if( t.message ) sb.append(t.message).append(' | ')
            t = t.cause
        }
        return sb.toString()
    }

    /**
     * Drive the script through a *real* {@link Session}, so the agent task is dispatched by the
     * real {@link nextflow.executor.local.AgentExecutor} into a real per-task work dir. The mock
     * harness cannot be used here: it substitutes every executor and never stages anything.
     */
    private Session runWithObserver(TraceObserverV2 probe, String text, Map config = [:]) {
        final workDir = Files.createTempDirectory('nxf-agent-path-io')
        def session = new Session([workDir: workDir.toString()] + config)
        session.setBinding(new ScriptBinding())
        session.init(null, null, null, null)
        // inject the probe into the private observersV2 list before ignition
        final f = Session.getDeclaredField('observersV2')
        f.setAccessible(true)
        final list = new ArrayList((List) f.get(session))
        list.add(probe)
        f.set(session, list)
        session.start()

        def loader = ScriptLoaderFactory.create(session)
        loader.parse(text)
        loader.runScript()

        session.fireDataflowNetwork()
        session.await()
        session.destroy()
        if( session.error )
            throw session.error
        return session
    }

}
