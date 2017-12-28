/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.executor
import java.nio.file.Files

import nextflow.Session
import nextflow.container.ContainerConfig
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondorExecutorTest extends Specification {


    def 'should return valid directives' () {
        given:
        def executor = [:] as CondorExecutor
        def task = new TaskRun()

        when:
        task.config = new TaskConfig()
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                queue
                                '''
                                .stripIndent().trim()


        when:
        task.config = new TaskConfig()
        task.config.cpus = 4
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                request_cpus = 4
                                machine_count = 1
                                queue
                                '''
                .stripIndent().trim()


        when:
        task.config = new TaskConfig()
        task.config.memory = '2 GB'
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                request_memory = 2 GB
                                queue
                                '''
                .stripIndent().trim()

        when:
        task.config = new TaskConfig()
        task.config.disk = '1 GB'
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                request_disk = 1 GB
                                queue
                                '''
                .stripIndent().trim()

        when:
        task.config = new TaskConfig()
        task.config.time = '1 h'
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                periodic_remove = (RemoteWallClockTime - CumulativeSuspensionTime) > 3600
                                queue
                                '''
                .stripIndent().trim()

        when:
        task.config = new TaskConfig()
        task.config.clusterOptions = 'alpha = 1; beta = 2'
        then:
        executor.getDirectives(task)
                .join('\n') ==  '''
                                universe = vanilla
                                executable = .command.run
                                log = .command.log
                                getenv = true
                                alpha = 1
                                beta = 2
                                queue
                                '''
                .stripIndent().trim()

    }


    def 'should return launch command line' () {

        given:
        def executor = [:] as CondorExecutor

        expect:
        executor.getSubmitCommandLine( Mock(TaskRun), null) == ['condor_submit', '--terse', '.command.condor']

    }

    def 'should parse job id' () {

        when:
        // executor stub object
        def executor = [:] as CondorExecutor
        then:
        executor.parseJobId( '15.0 - 15.0' ) == '15.0'

    }

    def 'should return the kill command' () {

        given:
        def executor = [:] as CondorExecutor

        expect:
        executor.getKillCommand() == ['condor_rm']

    }

    def 'should return kill command line' () {

        when:
        def executor = [:] as CondorExecutor
        then:
        executor.killTaskCommand('100.0') == ['condor_rm', '100.0']

    }

    def 'should parse queue output' () {

        given:
        def text = '''


        -- Submitter: 00979579dea7 : <10.0.0.17:9886?sock=58_3011_5> : 00979579dea7
         ID      OWNER            SUBMITTED     RUN_TIME ST PRI SIZE CMD
          20.0   submit          4/10 17:11   0+00:00:14 R  0   0.0  job.sh
          21.0   submit          4/10 17:11   0+00:00:00 I  0   0.0  job.sh
          22.0   submit          4/10 17:11   0+00:00:00 C  0   0.0  job.sh
          23.0   submit          4/10 17:11   0+00:00:00 E  0   0.0  job.sh
          24.0   submit          4/10 17:11   0+00:00:00 H  0   0.0  job.sh

        4 jobs; 0 completed, 0 removed, 3 idle, 1 running, 0 held, 0 suspended
        '''
        .stripIndent()

        when:
        def executor = [:] as CondorExecutor
        def status = executor.parseQueueStatus(text)

        then:
        status['20.0'] == AbstractGridExecutor.QueueStatus.RUNNING
        status['21.0'] == AbstractGridExecutor.QueueStatus.PENDING
        status['22.0'] == AbstractGridExecutor.QueueStatus.DONE
        status['23.0'] == AbstractGridExecutor.QueueStatus.ERROR
        status['24.0'] == AbstractGridExecutor.QueueStatus.HOLD
        status.size() == 5

    }

    def 'should parse empty status' () {

        given:
        def text = '''


        -- Submitter: a936e857590a : <172.17.0.2:57449> : a936e857590a
         ID      OWNER            SUBMITTED     RUN_TIME ST PRI SIZE CMD

        0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended
        '''
        .stripIndent()

        when:
        def executor = [:] as CondorExecutor
        def status = executor.parseQueueStatus(text)
        then:
        status.isEmpty()

        expect:
        executor.parseQueueStatus(null).isEmpty()
    }


    def 'should create condor submit file' () {

        given:
        def session = Mock(Session)
        session.getContainerConfig() >> new ContainerConfig(enabled:false)
        def folder = Files.createTempDirectory('test')
        def executor = [:] as CondorExecutor
        def task = new TaskRun(name: 'Hello', workDir: folder, script: 'echo Hello world!')
        task.config = Mock(TaskConfig)
        task.processor = Mock(TaskProcessor)
        task.processor.getProcessEnvironment() >> [:]
        task.processor.getConfig() >> Mock(ProcessConfig)
        task.processor.getSession() >> session

        /*
         * simple bash run
         */
        when:
        executor.createBashWrapperBuilder(task).build()

        then:
        folder.resolve('.command.run').exists()
        folder.resolve('.command.run').canExecute()

        folder.resolve('.command.sh').text == '''
                #!/bin/bash -ue
                echo Hello world!
                '''
                .stripIndent()
                .leftTrim()

        folder.resolve('.command.condor').text == '''
                    universe = vanilla
                    executable = .command.run
                    log = .command.log
                    getenv = true
                    queue
                    '''
                    .stripIndent()
                    .leftTrim()

    }

}
