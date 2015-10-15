/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.script

import java.nio.file.Files

import nextflow.Const
import nextflow.Session
import nextflow.scm.AssetManager
import nextflow.util.Duration
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
        def handlerInvoked = false
        def session = new Session([workflow: [onComplete: { -> handlerInvoked=true } ]  ])

        /*
         * script runner
         */
        def runner = Mock(ScriptRunner)
        runner.getScriptFile() >> script
        runner.fetchContainers() >> 'busybox/latest'
        runner.commandLine >> 'nextflow run -this -that'
        runner.session >> session

        when:
        def workflow = new WorkflowMetadata(runner)
        then:
        workflow.repository == 'https://github.com/nextflow-io/nextflow.git'
        workflow.commitId == commit.name().substring(0,10)
        workflow.revision == 'master'
        workflow.container == 'busybox/latest'
        workflow.localPath == dir
        workflow.startTime >= begin
        workflow.startTime <= new Date()
        workflow.endTime == null
        workflow.commandLine == 'nextflow run -this -that'
        workflow.nextflow.version == Const.APP_VER
        workflow.nextflow.build == Const.APP_BUILDNUM
        workflow.nextflow.timestamp == Const.APP_TIMESTAMP_UTC

        when:
        workflow.invokeOnComplete()
        then:
        workflow.endTime > workflow.startTime
        workflow.endTime <= new Date()
        workflow.duration == new Duration( workflow.endTime.time - workflow.startTime.time )
        handlerInvoked

        cleanup:
        dir?.deleteDir()
    }

}
