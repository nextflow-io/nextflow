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
import nextflow.trace.WorkflowStatsObserver
import nextflow.util.Duration
import nextflow.util.VersionNumber
import spock.lang.Specification
import spock.lang.Unroll
import test.TestHelper
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

    def 'should get a warning message' () {

        given:
        def session = new Session([process: ['$foo': [cpus:1], '$bar':[mem:'10GB']]])
        expect:
        session.validateConfig0(['foo','bar','baz']) == []
        session.validateConfig0(['foo','baz']) == ["There's no process matching config selector: bar -- Did you mean: baz?"]
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


    static Map cfg(String config) {
        new ConfigSlurper().parse(config).toMap()
    }


    def 'should fetch containers definition' () {

        String text

        when:
        text = '''
                process.container = 'beta'
                '''
        then:
        new Session(cfg(text)).fetchContainers() == 'beta'


        when:
        text = '''
                process {
                    withName:'proc1' { container = 'alpha' }
                    withName:'proc2' { container = 'beta' }
                }
                '''
        then:
        new Session(cfg(text)).fetchContainers() == ['proc1': 'alpha', 'proc2': 'beta']


        when:
        text = '''
                process {
                    withName:'proc1' { container = 'alpha' }
                    withName:'proc2' { container = 'beta' }
                }

                process.container = 'gamma'
                '''
        then:
        new Session(cfg(text)).fetchContainers() == ['proc1': 'alpha', 'proc2': 'beta', 'default': 'gamma']


        when:
        text = '''
                process.container = { "ngi/rnaseq:${workflow.getRevision() ?: 'latest'}" }
                '''

        def meta = Mock(WorkflowMetadata); meta.getRevision() >> '1.2'
        def session = new Session(cfg(text))
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
}
