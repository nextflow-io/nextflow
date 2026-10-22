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


package nextflow.trace

import java.nio.file.Files
import java.nio.file.Path

import groovyx.gpars.dataflow.Dataflow
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.operator.DataflowOperator
import nextflow.dag.DAG
import nextflow.file.FileHolder
import nextflow.processor.ErrorStrategy
import nextflow.processor.PublishDir
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.ProcessConfigV2
import nextflow.script.params.FileOutParam
import nextflow.trace.CleanupObserver.ProcessOutput
import nextflow.trace.CleanupObserver.ProcessState
import nextflow.trace.event.FilePublishEvent
import nextflow.trace.event.TaskEvent
import nextflow.trace.event.WorkflowOutputEvent
import spock.lang.Specification

class CleanupObserverTest extends Specification {

    DAG dag = new DAG()

    DAG.Vertex process(String name) {
        final vertex = new DAG.Vertex(dag, DAG.Type.PROCESS, name)
        vertex.process = Stub(TaskProcessor) { getName() >> name }
        dag.vertices << vertex
        return vertex
    }

    DAG.Vertex operator(List outputs) {
        final vertex = new DAG.Vertex(dag, DAG.Type.OPERATOR, 'map')
        vertex.operators = [ new DataflowOperator(Dataflow.DATA_FLOW_GROUP, [inputs: [new DataflowQueue()], outputs: outputs], { v -> }) ]
        dag.vertices << vertex
        return vertex
    }

    void edge(DAG.Vertex from, DAG.Vertex to, channel) {
        dag.edges << new DAG.Edge(dag, [from: from, to: to, channel: channel])
    }

    def 'should find direct and indirect consumers' () {
        given:
        def a1 = new DataflowQueue()
        def a2 = new DataflowQueue()
        def m1 = new DataflowQueue()
        def m2 = new DataflowQueue()
        def b1 = new DataflowQueue()
        def A = process('A')
        def B = process('B')
        def C = process('C')
        def op1 = operator([m1])
        def op2 = operator([m2])
        // A/1 -> map -> B
        edge(A, op1, a1)
        edge(op1, B, m1)
        // A/2 -> C, A/2 -> map -> (published)
        edge(A, C, a2)
        edge(A, op2, a2)
        edge(op2, null, m2)
        // B/1 -> (terminal)
        edge(B, null, b1)
        def outputs = [
            A: [ new ProcessOutput(a1, false), new ProcessOutput(a2, false) ],
            B: [ new ProcessOutput(b1, false) ],
            C: [],
        ]

        when:
        def states = CleanupObserver.getProcessStates(dag, outputs, [m2])

        then:
        states.A.consumers == [ ['B'] as Set, ['C'] as Set ]
        states.A.publishable == [ false, true ]
        states.B.consumers == [ [] as Set ]
        states.B.publishable == [ false ]
        states.C.consumers == []
    }

    def 'should mark consumers of forwarded inputs as indirect consumers' () {
        given:
        def a1 = new DataflowQueue()
        def b1 = new DataflowQueue()
        def A = process('A')
        def B = process('B')
        def C = process('C')
        // A/1 -> B, B/1 (forwards inputs) -> C, B/1 is published
        edge(A, B, a1)
        edge(B, C, b1)
        def outputs = [
            A: [ new ProcessOutput(a1, false) ],
            B: [ new ProcessOutput(b1, true) ],
            C: [],
        ]

        when:
        def states = CleanupObserver.getProcessStates(dag, outputs, [b1])

        then:
        states.A.consumers == [ ['B', 'C'] as Set ]
        states.A.publishable == [ true ]
        states.B.consumers == [ ['C'] as Set ]
        states.B.publishable == [ true ]
    }

    // -- event handling

    Path workDir

    def setup() {
        workDir = Files.createTempDirectory('cleanup')
    }

    def cleanup() {
        workDir.deleteDir()
    }

    Path file(String name) {
        final path = workDir.resolve(name)
        path.text = name
        return path
    }

    TaskRun task(String processName, List<Path> outputs, List<Path> inputs=[], boolean success=true, ErrorStrategy errorAction=null, List<PublishDir> publishDirs=[]) {
        final processor = Stub(TaskProcessor) { getName() >> processName }
        final param = Stub(FileOutParam) { getIndex() >> (short)0 }
        final inputFiles = inputs.collect { p -> new FileHolder(p) }
        Stub(TaskRun) {
            getName() >> processName
            getProcessor() >> processor
            getIndex() >> 1
            isSuccess() >> success
            isFailed() >> (errorAction != null)
            getErrorAction() >> errorAction
            getConfig() >> Stub(TaskConfig) { getPublishDir() >> publishDirs }
            getTargetDir() >> workDir
            getRestoreDir() >> null
            getOutputsByType(_) >> (outputs ? [(param): outputs] : [:])
            getInputFiles() >> inputFiles
            getInputFilesMap() >> inputFiles.collectEntries { h -> [h.stageName, h.storePath] }
        }
    }

    TaskEvent event(TaskRun task) {
        new TaskEvent(Stub(TaskHandler) { getTask() >> task })
    }

    CleanupObserver observer(Map<String,ProcessState> processes) {
        final observer = new CleanupObserver()
        observer.processes = processes
        return observer
    }

    def 'should delete a file once its consumers are done' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def B = task('B', [], [bam])
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when: 'producer completes'
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)

        when: 'consumer is pending'
        observer.onTaskPending(event(B))
        observer.onProcessClose(B.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)

        when: 'consumer completes'
        observer.onTaskComplete(event(B))
        observer.handleEvents()
        then:
        !Files.exists(bam)
        observer.paths.isEmpty()
        observer.consumerTasks.isEmpty()
    }

    def 'should keep deleting files when a file cannot be deleted' () {
        given:
        def readOnlyDir = Files.createDirectories(workDir.resolve('read-only'))
        def bam = readOnlyDir.resolve('a.bam')
        bam.text = 'a.bam'
        readOnlyDir.toFile().setWritable(false)
        def bai = file('a.bai')
        def A = task('A', [bam, bai])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ false ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.handleEvents()
        then:
        noExceptionThrown()
        Files.exists(bam)
        !Files.exists(bai)
        observer.paths.isEmpty()

        cleanup:
        readOnlyDir.toFile().setWritable(true)
    }

    def 'should not delete a file whose consumer was created before the producer was marked as completed' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def B = task('B', [], [bam])
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when: 'consumer is created and both processes are closed before the producer is marked as completed'
        observer.onTaskPending(event(B))
        observer.onProcessClose(A.processor)
        observer.onProcessClose(B.processor)
        observer.onTaskComplete(event(A))
        observer.handleEvents()
        then:
        Files.exists(bam)

        when: 'consumer completes'
        observer.onTaskComplete(event(B))
        observer.handleEvents()
        then:
        !Files.exists(bam)
        observer.consumerTasks.isEmpty()
    }

    def 'should not delete a file when the work directory path contains a symlink' () {
        given:
        def link = Files.createSymbolicLink(Files.createTempDirectory('cleanup-link').resolve('work'), workDir)
        def bam = link.resolve('a.bam')
        bam.text = 'a.bam'
        def A = task('A', [bam])
        def B = task('B', [], [bam])
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when: 'consumer is pending'
        observer.onTaskComplete(event(A))
        observer.onTaskPending(event(B))
        observer.onProcessClose(A.processor)
        observer.onProcessClose(B.processor)
        observer.handleEvents()
        then:
        B.inputFiles.first().storePath != bam
        Files.exists(bam)

        when: 'consumer completes'
        observer.onTaskComplete(event(B))
        observer.handleEvents()
        then:
        !Files.exists(bam)
        observer.consumerTasks.isEmpty()
    }

    def 'should not delete a file while its consumer process is open' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)

        when:
        observer.onProcessClose(Stub(TaskProcessor) { getName() >> 'B' })
        observer.handleEvents()
        then:
        !Files.exists(bam)
    }

    def 'should delete a publishable file only after it is published' () {
        given:
        def bam = file('a.bam')
        def target = workDir.resolve('published.bam')
        Files.createLink(target, bam)
        def A = task('A', [bam])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ true ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)

        when:
        observer.onFilePublish(new FilePublishEvent(bam, target))
        observer.handleEvents()
        then:
        !Files.exists(bam)
        Files.exists(target)
    }

    def 'should handle a file published before its task completes' () {
        given:
        def bam = file('a.bam')
        def target = workDir.resolve('published.bam')
        Files.createLink(target, bam)
        def other = file('collected.txt')
        def otherTarget = workDir.resolve('published.txt')
        Files.createLink(otherTarget, other)
        def A = task('A', [bam])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ true ]),
        ])

        when:
        observer.onFilePublish(new FilePublishEvent(bam, target))
        observer.onFilePublish(new FilePublishEvent(other, otherTarget))
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        !Files.exists(bam)
        observer.publishedOutputs == [other] as Set
    }

    def 'should delete an unpublished file once all workflow outputs are completed' () {
        given:
        def bam = file('a.bam')
        def bai = file('a.bai')
        def A = task('A', [bam])
        def B = task('B', [bai])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ true ]),
            B: new ProcessState([ [] as Set ], [ true ]),
        ])
        observer.pendingWorkflowOutputs = 2

        when: 'file is not published'
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onWorkflowOutput(new WorkflowOutputEvent('bam'))
        observer.handleEvents()
        then:
        Files.exists(bam)

        when: 'all workflow outputs are completed'
        observer.onWorkflowOutput(new WorkflowOutputEvent('summary'))
        observer.handleEvents()
        then:
        !Files.exists(bam)

        when: 'task completes after all workflow outputs are completed'
        observer.onTaskComplete(event(B))
        observer.onProcessClose(B.processor)
        observer.handleEvents()
        then:
        !Files.exists(bai)
    }

    def 'should not delete an unpublished file when workflow outputs are not tracked' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ true ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onWorkflowOutput(new WorkflowOutputEvent('bam'))
        observer.handleEvents()
        then:
        Files.exists(bam)
    }

    def 'should not delete a file that was published as a symlink' () {
        given:
        def bam = file('a.bam')
        def target = workDir.resolve('published.bam')
        Files.createSymbolicLink(target, bam)
        def A = task('A', [bam])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ true ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onFilePublish(new FilePublishEvent(bam, target))
        observer.handleEvents()
        then:
        Files.exists(bam)
        target.text == 'a.bam'
    }

    def 'should not delete the outputs of a failed task' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam], [], false)
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ false ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)
    }

    def 'should not delete the inputs of a failed task until it is retried' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def B1 = task('B', [], [bam], false, ErrorStrategy.RETRY)
        def B2 = task('B', [], [bam])
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when: 'consumer fails and will be retried'
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onTaskPending(event(B1))
        observer.onProcessClose(B1.processor)
        observer.onTaskComplete(event(B1))
        observer.handleEvents()
        then:
        Files.exists(bam)

        when: 'retry completes'
        observer.onTaskPending(event(B2))
        observer.onTaskComplete(event(B2))
        observer.handleEvents()
        then:
        !Files.exists(bam)
    }

    def 'should delete the inputs of a failed task that will not be retried' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam])
        def B = task('B', [], [bam], false, ErrorStrategy.IGNORE)
        def observer = observer([
            A: new ProcessState([ ['B'] as Set ], [ false ]),
            B: new ProcessState([], []),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onTaskPending(event(B))
        observer.onProcessClose(B.processor)
        observer.onTaskComplete(event(B))
        observer.handleEvents()
        then:
        !Files.exists(bam)
    }

    def 'should not track the outputs of a task that failed after a successful exit' () {
        given:
        def bam = file('a.bam')
        def A = task('A', [bam], [], true, ErrorStrategy.TERMINATE)
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ false ]),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.handleEvents()
        then:
        Files.exists(bam)
    }

    def 'should not delete a file published by publishDir until it is published' () {
        given:
        def bam = file('a.bam')
        def target = workDir.resolve('published.bam')
        Files.createLink(target, bam)
        def A = task('A', [bam], [], true, null, [PublishDir.create(path: workDir.resolve('results'), pattern: PATTERN)])
        def observer = observer([
            A: new ProcessState([ [] as Set ], [ false ]),
        ])
        observer.pendingWorkflowOutputs = 1

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onWorkflowOutput(new WorkflowOutputEvent('bam'))
        observer.handleEvents()
        then:
        Files.exists(bam) == PUBLISHED

        when:
        observer.onFilePublish(new FilePublishEvent(bam, target))
        observer.handleEvents()
        then:
        !Files.exists(bam)

        where:
        PATTERN     | PUBLISHED
        '*.bam'     | true
        '*.txt'     | false
    }

    def 'should mark consumers of an input file forwarded by a typed process as indirect consumers' () {
        given:
        def a = file('1.a')
        def d = file('1.d')
        def A = task('A', [a])
        def config = new ProcessConfigV2(Mock(BaseScript), 'D')
        config.getOutputs().addParam('out', List, null)
        def param = config.getOutputs().getParams()[0]
        def D = Stub(TaskRun) {
            getName() >> 'D'
            getProcessor() >> Stub(TaskProcessor) { getName() >> 'D'; getConfig() >> config }
            getIndex() >> 1
            isSuccess() >> true
            getConfig() >> Stub(TaskConfig) { getPublishDir() >> [] }
            getRestoreDir() >> null
            getOutputs() >> [(param): [a, d]]
            getOutputFiles() >> ([d] as Set)
            getInputFiles() >> [new FileHolder(a)]
        }
        def observer = observer([
            A: new ProcessState([ ['D'] as Set ], [ false ]),
            D: new ProcessState([ ['E'] as Set ], [ false ]),
            E: new ProcessState([], []),
        ])

        when:
        observer.onTaskComplete(event(A))
        observer.onProcessClose(A.processor)
        observer.onTaskPending(event(D))
        observer.onTaskComplete(event(D))
        observer.onProcessClose(D.processor)
        observer.handleEvents()
        then:
        Files.exists(a)
        Files.exists(d)

        when:
        observer.onProcessClose(Stub(TaskProcessor) { getName() >> 'E' })
        observer.handleEvents()
        then:
        !Files.exists(a)
        !Files.exists(d)
    }

}
