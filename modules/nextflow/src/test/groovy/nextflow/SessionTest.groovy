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

package nextflow

import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.attribute.PosixFilePermission
import java.util.concurrent.CountDownLatch
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit

import nextflow.config.Manifest
import nextflow.container.ContainerConfig
import nextflow.container.DockerConfig
import nextflow.container.PodmanConfig
import nextflow.container.SarusConfig
import nextflow.exception.AbortOperationException
import nextflow.file.FileHelper
import nextflow.script.ScriptFile
import nextflow.script.WorkflowMetadata
import nextflow.trace.TraceFileObserver
import nextflow.trace.TraceHelper
import nextflow.trace.TraceObserverV2
import nextflow.trace.WorkflowStatsObserver
import nextflow.util.Duration
import nextflow.util.VersionNumber
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper

import static test.ScriptHelper.*
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SessionTest extends Specification {


    def 'test baseDir and binDir'() {

        setup:
        def base = Files.createTempDirectory('test')
        def bin = base.resolve('bin'); bin.mkdir()

        when:
        def session = new Session()
        then:
        session.baseDir == null
        session.binDir == null

        when:
        session = new Session()
        session.baseDir = Paths.get('some/folder')
        then:
        session.baseDir == Paths.get('some/folder')
        session.binDir == null

        when:
        session = new Session()
        session.baseDir = base
        then:
        session.baseDir == base
        session.binDir.exists()

        cleanup:
        base.deleteDir()

    }

    def 'test add lib path'() {

        setup:
        def path1 = Files.createTempDirectory('test')
        def path2 = Files.createTempDirectory('test')

        when:
        def session = new Session()
        session.setLibDir( path1.toString() )
        then:
        session.getLibDir() == [ path1 ]


        when:
        session = new Session()
        session.setLibDir( "${path1}:${path2}" )
        then:
        session.getLibDir() == [ path1, path2 ]

        when:
        session = new Session()
        session.setBaseDir(Paths.get('/some/path'))
        then:
        session.getLibDir() == []

        when:
        def base = Files.createTempDirectory('test')
        base.resolve('lib').mkdir()
        session = new Session()
        session.setBaseDir(base)
        then:
        session.getLibDir() == [base.resolve('lib')]

        cleanup:
        base?.deleteDir()
        path1.deleteDir()
        path2.deleteDir()

    }

    def 'test cacheable property' () {
        when:
        def session = new Session()
        then:
        session.cacheable

        when:
        session = new Session([cacheable: false])
        then:
        !session.cacheable
    }

    def 'test create observers'() {
        given:
        TraceHelper.testTimestampFmt = '20221001'
        def session
        def result
        def observer

        when:
        session = [:] as Session
        result = session.createObserversV2()
        then:
        result.size()==1
        result.any { it instanceof WorkflowStatsObserver }

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, file:'name.txt']]
        result = session.createObserversV2()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == FileHelper.asPath('name.txt')
        observer.separator == '\t'

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, sep: 'x', fields: 'task_id,name,exit', file: 'alpha.txt']]
        result = session.createObserversV2()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == FileHelper.asPath('alpha.txt')
        observer.separator == 'x'
        observer.fields == ['task_id','name','exit']

        when:
        session = [:] as Session
        session.config = [trace: [sep: 'x', fields: 'task_id,name,exit']]
        result = session.createObserversV2()
        then:
        !result.any { it instanceof TraceFileObserver }

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, fields: 'task_id,name,exit,vmem']]
        result = session.createObserversV2()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == FileHelper.asPath('trace-20221001.txt')
        observer.separator == '\t'
        observer.fields == ['task_id','name','exit','vmem']

    }

    def 'should return absolute workDir' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def file = folder.resolve('pipeline.nf'); file.text = 'println "hello"'
        def script = new ScriptFile(file)

        when:
        def session = new Session([workDir: '../work'])
        session.init(script)

        then:
        session.binding != null
        session.baseDir == folder
        session.workDir.isAbsolute()
        !session.workDir.toString().contains('..')
        session.scriptName == 'pipeline.nf'
        session.classesDir.exists()
        session.observersV1 != null
        session.observersV2 != null
        session.workflowMetadata != null

        cleanup:
        session.classesDir?.deleteDir()

    }

    def 'should collect bin executable files' () {

        given:
        def folder = Files.createTempDirectory('test')
        Files.createFile(folder.resolve('foo.sh'))
        Files.createFile(folder.resolve('bar.sh'))
        Files.createFile(folder.resolve('baz.sh'))

        Files.setPosixFilePermissions(folder.resolve('foo.sh'), [PosixFilePermission.OWNER_READ, PosixFilePermission.OWNER_EXECUTE] as Set)
        Files.setPosixFilePermissions(folder.resolve('bar.sh'), [PosixFilePermission.OWNER_READ, PosixFilePermission.OWNER_EXECUTE] as Set)
        def session = [:] as Session

        when:
        def result = session.findBinEntries(folder)
        then:
        result.size() == 2
        result['foo.sh'] == folder.resolve('foo.sh')
        result['bar.sh'] == folder.resolve('bar.sh')

        cleanup:
        folder?.deleteDir()

    }

    def 'the agent orchestration pool is separate from the execution pool and grows on demand' () {
        given:
        def session = new Session([poolSize: 2])
        session.start()

        expect: 'the execution pool keeps its historical fixed shape, sized to poolSize'
        ((ThreadPoolExecutor) session.execService).getMaximumPoolSize() == 2

        and: 'orchestration draws from a DIFFERENT pool, so an agent cannot consume an execution thread'
        !session.getAgentExecService().is(session.execService)

        and: 'and holds no threads until an agent actually runs'
        ((ThreadPoolExecutor) session.getAgentExecService()).getPoolSize() == 0

        and: 'it is unbounded, so a blocked orchestrator can always be given a thread'
        ((ThreadPoolExecutor) session.getAgentExecService()).getMaximumPoolSize() == Integer.MAX_VALUE

        when: 'far more blocked orchestrators than the execution pool could ever admit'
        def pool = (ThreadPoolExecutor) session.getAgentExecService()
        int n = 50
        def started = new CountDownLatch(n)
        def release = new CountDownLatch(1)
        n.times { pool.submit({ started.countDown(); release.await() } as Runnable) }

        // sharing one fixed pool would admit only poolSize of these, and the tool sub-tasks they
        // block on would never get a thread -- the deadlock this partition removes
        then: 'all of them run concurrently'
        started.await(30, TimeUnit.SECONDS)
        pool.getPoolSize() >= n

        and: 'none of it consumed the execution pool'
        ((ThreadPoolExecutor) session.execService).getPoolSize() == 0

        cleanup:
        release.countDown()
        pool.shutdownNow()
        session.execService.shutdownNow()
    }

    def 'should get a warning message' () {

        given:
        def session = new Session([process: ['withName:foo': [cpus:1], 'withName:bar':[mem:'10GB']]])
        expect:
        session.validateConfig0(['foo','bar','baz']) == []
        session.validateConfig0(['foo','baz']) == ["There's no process matching config selector: bar -- Did you mean: baz?"]
    }

    def 'should validate agent selectors against the agent names' () {
        given:
        def session = new Session([agent: [
            model: 'openai/gpt-5-mini',
            'withName:critic': [cpus: 2],
            'withName:planer': [cpus: 4],
            'withLabel:reasoning': [cpus: 8] ]])

        expect: 'a matched agent selector is silent; the typo is reported as an agent, not a process'
        session.validateConfig0([], ['critic','planner']) == ["There's no agent matching config selector: planer -- Did you mean: planner?"]

        and: 'agent names never satisfy a `process` selector, and vice versa'
        new Session([process: ['withName:critic': [cpus:2]]]).validateConfig0([], ['critic'])
            == ["There's no process matching config selector: critic"]
        session.validateConfig0(['critic','planner'], []).size() == 2
    }

    @Unroll
    def 'should return engine type' () {
        given:
        def session =  new Session([(ENGINE): CONFIG])

        expect:
        session.containerConfig instanceof ContainerConfig
        session.containerConfig.enabled
        session.containerConfig.engine == ENGINE

        where:
        ENGINE         | CONFIG
        'docker'       | [enabled: true, x:'alpha', y: 'beta']
        'docker'       | [enabled: true, x:'alpha', y: 'beta', registry: 'd.reg']
        'podman'       | [enabled: true, x:'alpha', y: 'beta']
        'podman'       | [enabled: true, x:'alpha', y: 'beta', registry: 'd.reg']
        'sarus'        | [enabled: true, x:'delta', y: 'gamma']
        'shifter'      | [enabled: true, x:'delta', y: 'gamma']
        'singularity'  | [enabled: true, x:'delta', y: 'gamma']
        'charliecloud' | [enabled: true, x:'delta', y: 'gamma']
    }

    def 'should get config for specific engine' () {
        given:
        def config = [docker:[registry:'docker.io'], podman: [registry:'quay.io']]
        def session = new Session(config)

        expect:
        session.getContainerConfig(null) == new DockerConfig(registry:'docker.io')
        and:
        session.getContainerConfig('docker') == new DockerConfig(registry:'docker.io')
        and:
        session.getContainerConfig('podman') == new PodmanConfig(registry:'quay.io')
        and:
        session.getContainerConfig('sarus') == new SarusConfig([:])
    }

    @Unroll
    def 'should get config for conda environments' () {
        given:
        def session =  Spy(new Session([conda: CONFIG]))
        expect:
        session.condaConfig.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]

    }

    @Unroll
    def 'should get config for spack environments' () {
        given:
        def session =  Spy(new Session([spack: CONFIG]))
        expect:
        session.spackConfig.isEnabled() == EXPECTED

        where:
        EXPECTED    | CONFIG            | ENV
        false       | [:]               | [:]
        false       | [enabled: false]  | [:]
        true        | [enabled: true]   | [:]

    }

    def 'should get manifest object' () {

        given:
        def MAN = [author: 'pablo', nextflowVersion: '1.2.3', name: 'foo']

        when:
        def session = new Session([manifest: MAN])
        then:
        session.getManifest().with {
            author == 'pablo'
            nextflowVersion == '1.2.3'
            name == 'foo'
            description == null
        }
    }

    @Unroll
    def 'should check valid process name with selector=#SELECTOR' () {

        given:
        def session = new Session()

        when:
        def error = []
        session.checkValidProcessName(NAMES, SELECTOR, error)
        then:
        error[0] == MSG

        where:
        SELECTOR    | NAMES         | MSG
        'foo'       | ['foo','bar'] | null
        'bar'       | ['foo','bar'] | null
        'baz'       | ['foo','bar'] | "There's no process matching config selector: baz -- Did you mean: bar?"
        'ba.*'      | ['foo','bar'] | null
        'foo|bar'   | ['foo','bar'] | null
        'foo|bar'   | ['baz']       | "There's no process matching config selector: foo|bar"
    }

    def 'should fetch containers definition' () {

        String text

        when:
        text = '''
                process.container = 'beta'
                '''
        then:
        new Session(loadConfig(text)).fetchContainers() == 'beta'


        when:
        text = '''
                process {
                    withName:'proc1' { container = 'alpha' }
                    withName:'proc2' { container = 'beta' }
                }
                '''
        then:
        new Session(loadConfig(text)).fetchContainers() == ['proc1': 'alpha', 'proc2': 'beta']


        when:
        text = '''
                process {
                    withName:'proc1' { container = 'alpha' }
                    withName:'proc2' { container = 'beta' }
                }

                process.container = 'gamma'
                '''
        then:
        new Session(loadConfig(text)).fetchContainers() == ['proc1': 'alpha', 'proc2': 'beta', 'default': 'gamma']


        when:
        text = '''
                process.container = { "ngi/rnaseq:${workflow.getRevision() ?: 'latest'}" }
                '''

        def meta = Mock(WorkflowMetadata); meta.getRevision() >> '1.2'
        def session = new Session(loadConfig(text))
        session.binding.setVariable('workflow',meta)
        then:
        session.fetchContainers() == 'ngi/rnaseq:1.2'
    }

    @Unroll
    def 'should get module binaries status'() {
        given:
        def session = new Session()
        NextflowMeta.instance.moduleBinaries(MODE)

        expect:
        session.enableModuleBinaries() == EXPECTED

        cleanup:
        NextflowMeta.instance.moduleBinaries(false)

        where:
        MODE  | EXPECTED
        false | false
        true  | true

    }

    def 'should compute the auto resource labels lazily and memoize them' () {
        given:
        def session = new Session([tower: [autoLabels: 'runName']])
        // the Platform metadata is only filled in at `notifyFlowCreate` i.e. after the session is created
        def meta = Mock(WorkflowMetadata)
        session.@workflowMetadata = meta

        when:
        def labels = session.getAutoResourceLabels()
        then:
        // the metadata is only read on the first access
        (1.._) * meta.getRunName() >> 'crazy_darwin'
        and:
        labels == ['nextflow.io/runName': 'crazy_darwin']

        when:
        def again = session.getAutoResourceLabels()
        then:
        0 * meta.getRunName()
        and:
        again.is(labels)
    }

    def 'should return no auto resource labels when the option is not set' () {
        given:
        def session = new Session()
        session.@workflowMetadata = Mock(WorkflowMetadata)

        expect:
        session.getAutoResourceLabels() == [:]
    }

    @Unroll
    def 'should resolve the auto resource labels option scope' () {
        given:
        def session = new Session(CONFIG)
        session.@workflowMetadata = Mock(WorkflowMetadata) {
            getRunName() >> 'crazy_darwin'
            getProjectName() >> 'nf-core/rnaseq'
        }

        expect:
        session.getAutoResourceLabels().keySet() as List == EXPECTED

        where:
        CONFIG                                                                          | EXPECTED
        [tower: [autoLabels: 'runName']]                                                | ['nextflow.io/runName']
        [tower: [autoLabels: ['runName','projectName']]]                                | ['nextflow.io/projectName','nextflow.io/runName']
        [tower: [autoLabels: false]]                                                    | []
        [seqera: [executor: [autoLabels: 'runName']]]                                   | ['nextflow.io/runName']
        // the deprecated option wins when given, even as `false`
        [seqera: [executor: [autoLabels: 'projectName']], tower: [autoLabels:'runName']]| ['nextflow.io/projectName']
        [seqera: [executor: [autoLabels: false]], tower: [autoLabels: true]]            | []
        // ... and only when given
        [seqera: [executor: [endpoint: 'http://foo']], tower: [autoLabels: 'runName']]  | ['nextflow.io/runName']
    }

    def 'should report the offending option name for an invalid auto labels value' () {
        when:
        new Session([tower: [autoLabels: 'foo']]).getAutoResourceLabels()
        then:
        def e1 = thrown(IllegalArgumentException)
        e1.message.contains("'tower.autoLabels'")

        when:
        new Session([seqera: [executor: [autoLabels: 'foo']]]).getAutoResourceLabels()
        then:
        def e2 = thrown(IllegalArgumentException)
        e2.message.contains("'seqera.executor.autoLabels'")
    }

    def 'should notify flow complete only once when abort and destroy race' () {
        given:
        def observer = Mock(TraceObserverV2)
        def session = new Session()
        session.@observersV2 = [observer]

        when:
        // certain error conditions can abort the session on a separate thread
        // while the main thread winds it down via destroy()
        // both call shutdown0() -> notifyFlowComplete() concurrently
        def t1 = Thread.start { session.abort() }
        def t2 = Thread.start { session.destroy() }
        t1.join()
        t2.join()

        then:
        1 * observer.onFlowComplete()
    }
}
