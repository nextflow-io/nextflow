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

package nextflow.executor
import java.nio.file.Path
import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.MemoryUnit
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
                #BSUB -R "select[mem>=8192] rusage[mem=8192]"
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
                #BSUB -R "select[mem>=10] rusage[mem=10]"
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
                #BSUB -R "select[mem>=200] rusage[mem=200]"
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
                #BSUB -R "select[mem>=2048] rusage[mem=2048]"
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
                #BSUB -R "select[mem>=2048] rusage[mem=2048]"
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
                #BSUB -R "select[mem>=2048] rusage[mem=2048]"
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
                #BSUB -R "select[mem>=2048] rusage[mem=2048]"
                #BSUB -J nf-mapping_hola
                '''
                .stripIndent().leftTrim()

    }

    def testDiskResources() {
        given:
        def config = Mock(TaskConfig)
        def WORKDIR = Paths.get('/my/work')
        def lsf = [:] as LsfExecutor
        def task = Mock(TaskRun)

        when:
        def result = lsf.getDirectives(task)
        then:
        task.workDir >> WORKDIR
        task.config >> config
        task.name >> 'foo'
        config.getClusterOptionsAsList() >> []
        1 * config.getDisk() >> MemoryUnit.of('10GB')

        result.join(' ') == "-o $WORKDIR/.command.log -R select[tmp>=10240] rusage[tmp=10240] -J nf-foo"
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
                #BSUB -R "select[mem>=8192] rusage[mem=8192]"
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
                #BSUB -R "select[mem>=8192] rusage[mem=8192]"
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
                #BSUB -R "select[mem>=8192] rusage[mem=8192]"
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
                """\
                JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
                5085109 pluskal RUN   normal     it-c01b13   it-c01b11   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild/AutoBuild_run_1_/TEMP0/RUN_FILE_4 Oct 31 13:00
                5085585 pluskal RUN   normal     tak         it-c05b12   phenix.autobuild autobuild_config.eff Oct 31 14:05
                5085604 pluskal PEND  normal     it-c05b12   it-c05b12   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_1 Oct 31 14:44
                5085611 pluskal PSUSP normal     it-c05b12   it-c05b12   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_8 Oct 31 14:44
                5085107 pluskal EXIT  normal     it-c01b13   it-c01b07   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild/AutoBuild_run_1_/TEMP0/RUN_FILE_2 Oct 31 13:00
                5085607 pluskal UNKWN normal     it-c05b12   it-c01b07   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_4 Oct 31 14:44
                5085608 pluskal ZOMBI normal     it-c05b12   it-c01b05   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_5 Oct 31 14:44
                5085609 pluskal RUN   normal     it-c05b12   it-c01b14   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_6 Oct 31 14:44
                5085702 pluskal RUN   normal     tak         it-c01b14   ./assemble.nf Oct 31 15:08
                5085606 pluskal DONE  normal     it-c05b12   it-c01b12   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_3 Oct 31 14:44
                5085034 pluskal RUN   normal     tak         it-c01b13   phenix.autobuild autobuild_config.eff Oct 31 12:56
                5085612 pluskal RUN   normal     it-c05b12   it-c01b13   /lab/solexa_weng/xtal_data/Data_processing/Tomas_Pluskal/20181031_KLR1_autobuild_not_in_place/AutoBuild_run_1_/TEMP0/RUN_FILE_9 Oct 31 14:44
                5085703 pluskal RUN   normal     it-c01b14   it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31:it-r16u31 nf-trinityInchwormChrysalis Oct 31 15:08
                """.stripIndent().trim()


        when:
        def result = executor.parseQueueStatus(text)
        then:
        result['5085109'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085585'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085604'] == AbstractGridExecutor.QueueStatus.PENDING
        result['5085611'] == AbstractGridExecutor.QueueStatus.HOLD
        result['5085107'] == AbstractGridExecutor.QueueStatus.ERROR
        result['5085607'] == AbstractGridExecutor.QueueStatus.ERROR
        result['5085608'] == AbstractGridExecutor.QueueStatus.ERROR
        result['5085609'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085702'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085606'] == AbstractGridExecutor.QueueStatus.DONE
        result['5085034'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085612'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5085703'] == AbstractGridExecutor.QueueStatus.RUNNING
        result.size()==13
    }

    def 'should parse bjobs stats with extra headers' () {
        setup:
        def executor = [:] as LsfExecutor
        def TEXT = '''
            LSF is processing your request. Please wait ...
            LSF is processing your request. Please wait ...
            JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
            5157393 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_bd) Nov 14 13:00
            5157493 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_cn) Nov 14 13:00
            5157552 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_dl) Nov 14 13:00
            5157610 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_fe) Nov 14 13:00
            5157674 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_gf) Nov 14 13:00
            5157710 pluskal RUN   normal     it-c05b07   it-c05b10   *l_cmd_fv) Nov 14 13:00
            '''.stripIndent().trim()

        when:
        def result = executor.parseQueueStatus(TEXT)
        then:
        result['5157393'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5157493'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5157552'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5157610'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5157674'] == AbstractGridExecutor.QueueStatus.RUNNING
        result['5157710'] == AbstractGridExecutor.QueueStatus.RUNNING
        result.size() == 6 

    }


    def testQueueStatusCommand() {

        setup:
        def executor = [:] as LsfExecutor

        expect:
        executor.queueStatusCommand(null) == ['bjobs', '-w']
        executor.queueStatusCommand('long') == ['bjobs', '-w', '-q', 'long']

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

    def 'should parse col value' () {

        given:
        def HEADER = 'JOBID   USER    STAT  QUEUE'
        def executor = [:] as LsfExecutor
        def buf = new StringBuilder()

        expect:
        executor.col(HEADER,0,buf) == 'JOBID'
        executor.col(HEADER,8,buf) == 'USER'
        executor.col(HEADER,16,buf) == 'STAT'
    }
}
