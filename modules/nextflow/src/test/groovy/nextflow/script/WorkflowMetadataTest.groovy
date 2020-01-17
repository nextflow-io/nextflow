/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.nio.file.Files
import java.nio.file.Paths
import java.time.OffsetDateTime

import nextflow.Const
import nextflow.Session
import nextflow.scm.AssetManager
import nextflow.trace.TraceRecord
import nextflow.trace.WorkflowStats
import nextflow.trace.WorkflowStatsObserver
import nextflow.util.Duration
import nextflow.util.VersionNumber
import org.eclipse.jgit.api.Git
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowMetadataTest extends Specification {

    def 'should populate workflow object' () {

        given:
        final begin = OffsetDateTime.now()
        def dir = Files.createTempDirectory('test')
        /*
         * create the github repository
         */
        dir.resolve('main.nf').text = "println 'Hello world'"
        dir.resolve('nextflow.config').text = 'manifest {  }'

        def init = Git.init()
        def repo = init.setDirectory( dir.toFile() ).call()
        repo.add().addFilepattern('.').call()
        def commit = repo.commit().setAll(true).setMessage('First commit').call()
        repo.close()

        // append fake remote data
        dir.resolve('.git/config') << '''
            [remote "origin"]
                url = https://github.com/nextflow-io/rnaseq-nf.git
                fetch = +refs/heads/*:refs/remotes/origin/*
            [branch "master"]
                remote = origin
                merge = refs/heads/master
            '''
                .stripIndent()

        /*
         * create ScriptFile object
         */
        def manager = new AssetManager()
                    .setLocalPath(dir.toFile())
                    .setProject('nextflow-io/rnaseq-nf')
        def script = manager.getScriptFile()

        /*
         * config file onComplete handler
         */
        def handlerInvoked
        def config = [workflow: [onComplete: { -> handlerInvoked=workflow.commandLine } ],
                    docker:[enabled:true],
                    manifest: [version: '1.0.0', nextflowVersion: '>=0.31.1']]
        Session session = Spy(Session, constructorArgs: [config])
        session.configFiles >> [Paths.get('foo'), Paths.get('bar')]
        session.getStatsObserver() >> Mock(WorkflowStatsObserver) { getStats() >> new WorkflowStats() }
        session.fetchContainers() >> 'busybox/latest'
        session.commandLine >> 'nextflow run -this -that'

        when:
        def metadata = new WorkflowMetadata(session, script)
        session.binding.setVariable('workflow',metadata)
        then:
        metadata.scriptId == '0e44b16bdeb8ef9f4d8aadcf520e717d'
        metadata.scriptFile == manager.scriptFile.main
        metadata.scriptName == 'main.nf'
        metadata.repository == 'https://github.com/nextflow-io/rnaseq-nf.git'
        metadata.commitId == commit.name()
        metadata.revision == 'master'
        metadata.container == 'busybox/latest'
        metadata.projectDir == dir
        metadata.projectName == 'nextflow-io/rnaseq-nf'
        metadata.start >= begin
        metadata.start <= OffsetDateTime.now()
        metadata.complete == null
        metadata.commandLine == 'nextflow run -this -that'
        metadata.nextflow.version == new VersionNumber(Const.APP_VER)
        metadata.nextflow.build == Const.APP_BUILDNUM
        metadata.nextflow.timestamp == Const.APP_TIMESTAMP_UTC
        metadata.profile == 'standard'
        metadata.sessionId == session.uniqueId
        metadata.runName == session.runName
        metadata.containerEngine == 'docker'
        metadata.configFiles == [Paths.get('foo').toAbsolutePath(), Paths.get('bar').toAbsolutePath()]
        metadata.resume == false
        metadata.userName == System.getProperty('user.name')
        metadata.homeDir == Paths.get(System.getProperty('user.home'))
        metadata.manifest.version == '1.0.0'
        metadata.manifest.nextflowVersion == '>=0.31.1'

        when:
        metadata.invokeOnComplete()
        then:
        metadata.complete > metadata.start
        metadata.complete <= OffsetDateTime.now()
        metadata.duration == Duration.between(metadata.start, metadata.complete)
        handlerInvoked == metadata.commandLine


        when:
        session.profile >> 'foo_profile'
        metadata = new WorkflowMetadata(session, script)
        then:
        metadata.profile == 'foo_profile'

        cleanup:
        dir?.deleteDir()
    }

    def foo_test_method() {
        return 'foo_value'
    }

    def 'should access workflow script variables onComplete' () {

        given:
        def file = TestHelper.createInMemTempFile('main.nf', "println 'Hello world'")
        def script = new ScriptFile(file)

        def session = Spy(Session)
        session.getStatsObserver() >> Mock(WorkflowStatsObserver) { getStats() >> new WorkflowStats() }
        
        def metadata = new WorkflowMetadata(session, script)

        session.binding.setVariable('value_a', 1)
        session.binding.setVariable('value_b', 2)
        session.binding.setVariable('workflow', metadata)
        session.binding.setParams(foo: 'Hello', bar: 'world')

        def result1
        def result2
        def result3
        def result4
        def result5
        def result6
        def result7

        def handler = {
            result1 = workflow.commandLine   // workflow property
            result2 = workflow      // workflow object itself
            result3 = value_a       // variable in the session binding
            result4 = events        // workflow private field, should not be accessed
            result5 = xyz           // unknown field, should return null
            result6 = foo_test_method()
            result7 = "$params.foo $params.bar"
        }

        when:
        metadata.onComplete(handler)
        metadata.invokeOnComplete()

        then:
        result1 == metadata.commandLine
        result2 == metadata
        result3 == 1
        result4 == null
        result5 == null
        result6 == 'foo_value'
        result7 == 'Hello world'


    }

    def 'should access workflow script variables onError' () {

        given:
        def file = TestHelper.createInMemTempFile('main.nf', "println 'Hello world'")
        def script = new ScriptFile(file)

        def session = Spy(Session)
        def metadata = new WorkflowMetadata(session, script)
        
        session.binding.setVariable('value_a', 1)
        session.binding.setVariable('value_b', 2)
        session.binding.setVariable('workflow', metadata)
        session.binding.setParams(foo: 'Hello', bar: 'world')
        
        session.getStatsObserver() >> Mock(WorkflowStatsObserver) { getStats() >> new WorkflowStats() }

        def result1
        def result2
        def result3
        def result4
        def result5
        def result6
        def result7
        def result8

        def handler = {
            result1 = workflow.commandLine   // workflow property
            result2 = workflow      // workflow object itself
            result3 = value_a       // variable in the session binding
            result4 = events        // workflow private field, should not be accessed
            result5 = xyz           // unknown field, should return null
            result6 = foo_test_method()
            result7 = "$params.foo $params.bar"
            result8 = workflow.success
        }

        when:
        metadata.onError(handler)
        metadata.invokeOnError(Mock(TraceRecord))

        then:
        result1 == metadata.commandLine
        result2 == metadata
        result3 == 1
        result4 == null
        result5 == null
        result6 == 'foo_value'
        result7 == 'Hello world'
        result8 == false

    }

    def 'should convert to map object' () {

        given:
        def script = Mock(ScriptFile)
        script.getScriptId() >> '123'
        script.getCommitId() >> 'abcd'
        
        def session = Spy(Session)
        def metadata = new WorkflowMetadata(session, script)
        and:
        session.getStatsObserver() >> Mock(WorkflowStatsObserver) { getStats() >> new WorkflowStats() }

        when:
        def result = metadata.toMap()
        println result
        then:
        result.scriptId == '123'
        result.commitId == 'abcd'
        result.profile == 'standard'
    }

}
