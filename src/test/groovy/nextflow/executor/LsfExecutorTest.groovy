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
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LsfExecutorTest extends Specification {

    def testCommandLine() {

        when:
        def executor = [:] as LsfExecutor
        then:
        executor.getSubmitCommandLine(Mock(TaskRun), null) == ['bsub']

    }


    def testHeaders() {

        setup:
        // LSF executor
        def executor = [:] as LsfExecutor
        executor.session = new Session()

        // mock process
        def proc = Mock(TaskProcessor)
        // process name
        proc.getName() >> 'task'

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/scratch')
        task.name = 'mapping hola'

        when:
        task.config = new TaskConfig()
        // config
        task.config.queue = 'bsc_ls'
        task.config.clusterOptions = "-x 1 -R \"span[ptile=2]\""
        task.config.cpus = '2'
        task.config.time = '1h 30min'
        task.config.memory = '8GB'

        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q bsc_ls
                #BSUB -n 2
                #BSUB -R "span[hosts=1]"
                #BSUB -W 01:30
                #BSUB -M 4096
                #BSUB -R "rusage[mem=8192]"
                #BSUB -J nf-mapping_hola
                #BSUB -x 1
                #BSUB -R "span[ptile=2]"
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.cpus = 1
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q alpha
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'alpha'
        task.config.cpus = 1
        task.config.time = '1min'
        task.config.memory = '10MB'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q alpha
                #BSUB -W 00:01
                #BSUB -M 10
                #BSUB -R "rusage[mem=10]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'gamma'
        task.config.cpus = 1
        task.config.time = '4h'
        task.config.memory = '200MB'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q gamma
                #BSUB -W 04:00
                #BSUB -M 200
                #BSUB -R "rusage[mem=200]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'gamma'
        task.config.cpus = 4
        task.config.memory = '2GB'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q gamma
                #BSUB -n 4
                #BSUB -R "span[hosts=1]"
                #BSUB -M 512
                #BSUB -R "rusage[mem=2048]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()

        when:
        task.config = new TaskConfig()
        task.config.queue = 'gamma'
        task.config.cpus = 4
        task.config.memory = '2GB'
        task.config.time = '1d'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q gamma
                #BSUB -n 4
                #BSUB -R "span[hosts=1]"
                #BSUB -W 24:00
                #BSUB -M 512
                #BSUB -R "rusage[mem=2048]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'gamma'
        task.config.cpus = 8
        task.config.memory = '2GB'
        task.config.time = '2d'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q gamma
                #BSUB -n 8
                #BSUB -R "span[hosts=1]"
                #BSUB -W 48:00
                #BSUB -M 256
                #BSUB -R "rusage[mem=2048]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()


        when:
        task.config = new TaskConfig()
        task.config.queue = 'delta'
        task.config.memory = '2GB'
        task.config.time = '2d 12h 5m'
        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q delta
                #BSUB -W 60:05
                #BSUB -M 2048
                #BSUB -R "rusage[mem=2048]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()

    }

    def testPerJobMemLimit() {
        setup:
        // mock process
        def proc = Mock(TaskProcessor)
        // process name
        proc.getName() >> 'task'

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/scratch')
        task.name = 'mapping hola'

        // config
        task.config = new TaskConfig()
        task.config.queue = 'bsc_ls'
        task.config.cpus = 4
        task.config.memory = '8GB'

        when:
        // LSF executor
        def executor = [:] as LsfExecutor
        executor.session = new Session()

        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q bsc_ls
                #BSUB -n 4
                #BSUB -R "span[hosts=1]"
                #BSUB -M 2048
                #BSUB -R "rusage[mem=8192]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()

        when:
        // LSF executor
        executor = [:] as LsfExecutor
        executor.session = new Session([executor: [perJobMemLimit: true]])

        then:
        executor.getHeaders(task) == '''
                #BSUB -o /scratch/.command.log
                #BSUB -q bsc_ls
                #BSUB -n 4
                #BSUB -R "span[hosts=1]"
                #BSUB -M 8192
                #BSUB -R "rusage[mem=8192]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()

    }

    def testWorkDirWithBlanks() {

        setup:
        // LSF executor
        def executor = Spy(LsfExecutor)
        executor.session = new Session()

        // mock process
        def proc = Mock(TaskProcessor)
        // process name
        proc.getName() >> 'task'

        // task object
        def task = new TaskRun()
        task.processor = proc
        task.workDir = Paths.get('/scratch/some data/path')
        task.name = 'mapping hola'

        when:
        task.config = new TaskConfig()
        // config
        task.config.queue = 'bsc_ls'
        task.config.clusterOptions = "-x 1 -R \"span[ptile=2]\""
        task.config.cpus = '2'
        task.config.time = '1h 30min'
        task.config.memory = '8GB'


        then:
        executor.getHeaders(task) == '''
                #BSUB -o "/scratch/some data/path/.command.log"
                #BSUB -q bsc_ls
                #BSUB -n 2
                #BSUB -R "span[hosts=1]"
                #BSUB -W 01:30
                #BSUB -M 4096
                #BSUB -R "rusage[mem=8192]"
                #BSUB -J nf-mapping_hola
                #BSUB -x 1
                #BSUB -R "span[ptile=2]"
                '''
                .stripIndent().leftTrim()


    }

    def testParseJobId() {

        when:
        // executor stub object
        def executor = [:] as LsfExecutor
        then:
        executor.parseJobId( 'Job <2329803> is submitted to default queue <research-rh6>.' ) == '2329803'

    }

    def testKillCommand() {
        when:
        // executor stub object
        def executor = [:] as LsfExecutor
        then:
        executor.killTaskCommand('12345').join(' ') == 'bkill 12345'

    }

    def testQstatCommand() {

        setup:
        def executor = [:] as LsfExecutor
        def text =
                """
                6795348,RUN,Feb 17 13:26
                6795349,RUN,Feb 17 13:26
                6795351,PEND,Feb 17 13:26
                6795353,PSUSP,Feb 17 13:26
                6795354,EXIT,Feb 17 13:26
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result.size() == 5
        result['6795348'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['6795349'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['6795351'] == AbstractGridExecutor.QueueStatus.PENDING
        result['6795353'] == AbstractGridExecutor.QueueStatus.HOLD
        result['6795354'] == AbstractGridExecutor.QueueStatus.ERROR

    }


    def testQueueStatusCommand() {

        setup:
        def executor = [:] as LsfExecutor

        expect:
        executor.queueStatusCommand(null) == ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader']
        executor.queueStatusCommand('long') == ['bjobs', '-o',  'JOBID STAT SUBMIT_TIME delimiter=\',\'', '-noheader', '-q', 'long']

    }



    def testWrapString() {

        given:
        def executor = [:] as LsfExecutor

        expect:
        executor.wrapHeader('') == ''
        executor.wrapHeader('hello') == 'hello'
        executor.wrapHeader('span[ptile=number]') == '"span[ptile=number]"'
        executor.wrapHeader('hello world') == '"hello world"'

        executor.pipeLauncherScript() == true
        executor.getSubmitCommandLine(Mock(TaskRun), Mock(Path)) == ['bsub']
    }
}
