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

package nextflow.script

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.Const
import nextflow.Session
import nextflow.scm.AssetManager
import nextflow.trace.TraceRecord
import nextflow.util.Duration
import nextflow.util.VersionNumber
import org.eclipse.jgit.api.Git
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WorkflowMetadataTest extends Specification {

    def 'should populate workflow object' () {

        given:
        final begin = new Date()
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
                url = https://github.com/nextflow-io/nextflow.git
                fetch = +refs/heads/*:refs/remotes/origin/*
            [branch "master"]
                remote = origin
                merge = refs/heads/master
            '''
                .stripIndent()

        /*
         * create ScriptFile object
         */
        def manager = new AssetManager().setLocalPath(dir.toFile())
        def script = manager.getScriptFile()

        /*
         * config file onComplete handler
         */
        def handlerInvoked
        def session = new Session([workflow: [onComplete: { -> handlerInvoked=workflow.commandLine } ],
                                   docker:[enabled:true],
                                   manifest: [version: '1.0.0', nextflowVersion: '>=0.31.1']])
        session.configFiles = [Paths.get('foo'), Paths.get('bar')]
        /*
         * script runner
         */
        def runner = Mock(ScriptRunner)
        runner.getScriptFile() >> script
        runner.fetchContainers() >> 'busybox/latest'
        runner.commandLine >> 'nextflow run -this -that'
        runner.session >> session

        when:
        def metadata = new WorkflowMetadata(runner)
        session.binding.setVariable('workflow',metadata)
        then:
        metadata.scriptId == '0e44b16bdeb8ef9f4d8aadcf520e717d'
        metadata.scriptFile == manager.scriptFile.main
        metadata.scriptName == 'main.nf'
        metadata.repository == 'https://github.com/nextflow-io/nextflow.git'
        metadata.commitId == commit.name()
        metadata.revision == 'master'
        metadata.container == 'busybox/latest'
        metadata.projectDir == dir
        metadata.start >= begin
        metadata.start <= new Date()
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
        metadata.complete <= new Date()
        metadata.duration == new Duration( metadata.complete.time - metadata.start.time )
        handlerInvoked == metadata.commandLine


        when:
        runner.profile >> 'foo_profile'
        metadata = new WorkflowMetadata(runner)
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
        def dir = Files.createTempDirectory('test')
        dir.resolve('main.nf').text = "println 'Hello world'"
        def script = new ScriptFile(dir.resolve('main.nf').toFile())

        def session = new Session()

        def runner = Mock(ScriptRunner)
        runner.getScriptFile() >> script
        runner.fetchContainers() >> 'busybox/latest'
        runner.commandLine >> 'nextflow run -this -that'
        runner.session >> session

        when:
        def metadata = new WorkflowMetadata(runner)
        session.binding.setVariable('value_a', 1)
        session.binding.setVariable('value_b', 2)
        session.binding.setVariable('workflow', metadata)

        def result1
        def result2
        def result3
        def result4
        def result5
        def result6

        def handler = {
            result1 = workflow.commandLine   // workflow property
            result2 = workflow      // workflow object itself
            result3 = value_a       // variable in the session binding
            result4 = events        // workflow private field, should not be accessed
            result5 = xyz           // unknown field, should return null
            result6 = foo_test_method()
        }
        metadata.onComplete(handler)
        metadata.invokeOnComplete()

        then:
        result1 == metadata.commandLine
        result2 == metadata
        result3 == 1
        result4 == null
        result5 == null
        result6 == 'foo_value'

        cleanup:
        dir?.deleteDir()

    }

    def 'should access workflow script variables onError' () {

        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('main.nf').text = "println 'Hello world'"
        def script = new ScriptFile(dir.resolve('main.nf').toFile())

        def session = new Session()

        def runner = Mock(ScriptRunner)
        runner.getScriptFile() >> script
        runner.fetchContainers() >> 'busybox/latest'
        runner.commandLine >> 'nextflow run -this -that'
        runner.session >> session

        when:
        def metadata = new WorkflowMetadata(runner)
        session.binding.setVariable('value_a', 1)
        session.binding.setVariable('value_b', 2)
        session.binding.setVariable('workflow', metadata)

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
            result7 = workflow.success
        }
        metadata.onError(handler)
        metadata.invokeOnError(Mock(TraceRecord))

        then:
        result1 == metadata.commandLine
        result2 == metadata
        result3 == 1
        result4 == null
        result5 == null
        result6 == 'foo_value'
        result7 == false

        cleanup:
        dir?.deleteDir()

    }


}
