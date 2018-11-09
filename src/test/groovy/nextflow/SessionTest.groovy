/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import nextflow.container.ContainerConfig
import nextflow.trace.StatsObserver
import nextflow.trace.TraceFileObserver
import nextflow.util.Duration
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


    def 'test get queue size'() {

        when:
        def session = [:] as Session
        session.config = [ executor:['$sge':[queueSize: 123] ] ]
        then:
        session.getQueueSize('sge', 1) == 123
        session.getQueueSize('xxx', 1) == 1
        session.getQueueSize(null, 1) == 1

        when:
        def session2 = [:] as Session
        session2.config = [ executor:[ queueSize: 321, '$sge':[queueSize:789] ] ]
        then:
        session2.getQueueSize('sge', 2) == 789
        session2.getQueueSize('xxx', 2) == 321
        session2.getQueueSize(null, 2) == 321


        when:
        def session3 = [:] as Session
        session3.config = [ executor: 'sge' ]
        then:
        session3.getQueueSize('sge', 1) == 1
        session3.getQueueSize('xxx', 2) == 2
        session3.getQueueSize(null, 3) == 3


    }

    def 'test get poll interval'() {

        when:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[pollInterval: 345] ] ]
        then:
        session1.getPollInterval('sge').toMillis() == 345
        session1.getPollInterval('xxx').toMillis() == 1_000
        session1.getPollInterval(null).toMillis() == 1_000
        session1.getPollInterval(null, 2_000 as Duration).toMillis() == 2_000

        when:
        def session2 = [:] as Session
        session2.config = [ executor:[ pollInterval: 321, '$sge':[pollInterval:789] ] ]
        then:
        session2.getPollInterval('sge').toMillis() == 789
        session2.getPollInterval('xxx').toMillis() == 321
        session2.getPollInterval(null).toMillis() == 321

        when:
        def session3 = [:] as Session
        session3.config = [ executor: 'lsf' ]
        then:
        session3.getPollInterval('sge', 33 as Duration ).toMillis() == 33
        session3.getPollInterval('xxx', 44 as Duration ).toMillis() == 44
        session3.getPollInterval(null, 55 as Duration).toMillis() == 55

    }

    def 'test get exit read timeout'() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[exitReadTimeout: '5s'] ] ]

        expect:
        session1.getExitReadTimeout('sge') == '5sec' as Duration
        session1.getExitReadTimeout('lsf', '3sec' as Duration) == '3sec' as Duration

    }

    def 'test get queue stat interval'() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[queueStatInterval: '4sec'] ] ]

        expect:
        session1.getQueueStatInterval('sge') == '4sec' as Duration
        session1.getQueueStatInterval('lsf', '1sec' as Duration) == '1sec' as Duration

    }

    def 'test monitor dump interval'() {

        setup:
        def session1 = [:] as Session
        session1.config = [ executor:['$sge':[dumpInterval: '6sec'] ] ]

        expect:
        session1.getMonitorDumpInterval('sge') == '6sec' as Duration
        session1.getMonitorDumpInterval('lsf', '2sec' as Duration) == '2sec' as Duration

    }

    def 'test get exec config prop'() {

        when:
        def session = [:] as Session
        session.config = [ executor: [x:123, y:222, '$hazelcast': [y:333] ] ]
        then:
        session.getExecConfigProp( 'hazelcast', 'x', null ) == 123
        session.getExecConfigProp( 'hazelcast', 'y', null ) == 333
        session.getExecConfigProp( 'local', 'y', null ) == 222
        session.getExecConfigProp( 'local', 'y', 'beta') == 222
        session.getExecConfigProp( 'hazelcast', 'z', null ) ==  null
        session.getExecConfigProp( 'hazelcast', 'z', 'alpha') == 'alpha'
        session.getExecConfigProp( 'hazelcast', 'z', 'alpha', [NXF_EXECUTOR_Z:'hola']) == 'hola'
        session.getExecConfigProp( 'hazelcast', 'p.q.z', null, [NXF_EXECUTOR_P_Q_Z:'hello']) == 'hello'
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

    def 'test create observers'() {

        def session
        def result
        def observer

        when:
        session = [:] as Session
        result = session.createObservers()
        then:
        result.size()==1
        result.any { it instanceof StatsObserver }

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, file:'name.txt']]
        result = session.createObservers()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == Paths.get('name.txt').complete()
        observer.separator == '\t'

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, sep: 'x', fields: 'task_id,name,exit', file: 'alpha.txt']]
        result = session.createObservers()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == Paths.get('alpha.txt').complete()
        observer.separator == 'x'
        observer.fields == ['task_id','name','exit']

        when:
        session = [:] as Session
        session.config = [trace: [sep: 'x', fields: 'task_id,name,exit']]
        result = session.createObservers()
        then:
        !result.any { it instanceof TraceFileObserver }

        when:
        session = [:] as Session
        session.config = [trace: [enabled: true, fields: 'task_id,name,exit,vmem']]
        result = session.createObservers()
        observer = result[1] as TraceFileObserver
        then:
        result.size() == 2
        observer.tracePath == Paths.get('trace.txt').complete()
        observer.separator == '\t'
        observer.fields == ['task_id','name','exit','vmem']


    }

    def 'should return absolute workDir' () {

        given:
        def folder = TestHelper.createInMemTempDir()
        def script = folder.resolve('pipeline.nf')

        when:
        def session = new Session([workDir: '../work'])
        session.init(script)
        then:
        session.baseDir == folder
        session.workDir.isAbsolute()
        !session.workDir.toString().contains('..')

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
        session.validateConfig0(['foo','baz']) == ['The config file defines settings for an unknown process: bar -- Did you mean: baz?']
    }

    @Unroll
    def 'should return engine type' () {
        given:
        def session =  new Session([(engine): config])

        expect:
        session.containerConfig == config as ContainerConfig
        session.containerConfig.enabled
        session.containerConfig.engine == engine

        where:
        engine      | config
        'docker'    | [enabled: true, x:'alpha', y: 'beta']
        'docker'    | [enabled: true, x:'alpha', y: 'beta', registry: 'd.reg']
        'udocker'   | [enabled: true, x:'alpha', y: 'beta']
        'shifter'   | [enabled: true, x:'delta', y: 'gamma']
    }

    def 'should get manifest object' () {

        given:
        def MAN = [author: 'pablo', nextflowVersion: '1.2.3', name: 'foo']

        when:
        def session = new Session([manifest: MAN])
        then:
        session.manifest.with {
            author == 'pablo'
            nextflowVersion == '1.2.3'
            name == 'foo'
            description == null
        }
    }

    def 'should get config attribute' () {

        given:
        def session = Spy(Session)

        when:
        def result = session.getConfigAttribute('alpha', 'hello')
        then:
        result == 'hello'

        when:
        result = session.getConfigAttribute('delta', 'hello')
        then:
        session.getConfig() >> [delta: '1234']
        result == '1234'

        when:
        result = session.getConfigAttribute('omega', 'hello')
        then:
        session.getSystemEnv() >> [NXF_OMEGA: '6789']
        result == '6789'
    }

    def 'should get config nested attribute' () {

        given:
        def session = Spy(Session)

        when:
        def result = session.getConfigAttribute('alpha.beta.delta', 'hello')
        then:
        result == 'hello'

        when:
        result = session.getConfigAttribute('alpha.beta.gamma', 'hello')
        then:
        session.getConfig() >> [alpha: [beta: [gamma: 'abc']]]
        result == 'abc'

        when:
        result = session.getConfigAttribute('alpha.beta.omega', 'hello')
        then:
        session.getSystemEnv() >> [NXF_ALPHA_BETA_OMEGA: 'OK']
        result == 'OK'

    }


}
