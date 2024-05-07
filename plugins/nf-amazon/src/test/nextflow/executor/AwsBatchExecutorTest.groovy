/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.executor

import java.nio.file.Path

import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.aws.batch.AwsBatchExecutor
import nextflow.cloud.aws.batch.AwsOptions
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchExecutorTest extends Specification {

    def 'should check is fusion' () {
        given:
        SysEnv.push(ENV)
        and:
        def sess = Mock(Session) {
            getConfig() >> CONFIG
        }
        def executor = new AwsBatchExecutor(session: sess)

        expect:
        executor.isFusionEnabled() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        CONFIG                      | ENV                           | EXPECTED
        [:]                         | [:]                           | false
        [fusion:[enabled: true]]    | [:]                           | true
        [fusion:[enabled: false]]   | [FUSION_ENABLED:'true']   | false     // <-- config has priority
        [:]                         | [FUSION_ENABLED:'true']   | true
        [:]                         | [FUSION_ENABLED:'false']  | false

    }

    def 'should validate shouldDeleteJob method' () {
        given:
        def executor = Spy(AwsBatchExecutor)

        expect:
        executor.shouldDeleteJob('job-1')
        executor.shouldDeleteJob('job-2')
        executor.shouldDeleteJob('job-3')
        and:
        !executor.shouldDeleteJob('job-1')
        !executor.shouldDeleteJob('job-1')
        !executor.shouldDeleteJob('job-2')
        !executor.shouldDeleteJob('job-2')
        !executor.shouldDeleteJob('job-3')
        !executor.shouldDeleteJob('job-3')
    }

    def 'should get array index variable and start' () {
        given:
        def executor = Spy(AwsBatchExecutor)
        expect:
        executor.getArrayIndexName() == 'AWS_BATCH_JOB_ARRAY_INDEX'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(AwsBatchExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:
        JOB_ID      | TASK_INDEX    | EXPECTED
        'foo'       | 1             | 'foo:1'
        'bar'       | 2             | 'bar:2'
    }

    protected Path s3(String path) {
        S3PathFactory.parse('s3:/' + path)
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(AwsBatchExecutor) {
            isFusionEnabled()>>FUSION
            isWorkDirDefaultFS()>>DEFAULT_FS
        }
        and:
        def handler = Mock(TaskHandler) {
            getTask() >> Mock(TaskRun) { getWorkDir() >> WORK_DIR }
        }
        expect:
        executor.getArrayWorkDir(handler) == EXPECTED

        where:
        FUSION  | DEFAULT_FS  | WORK_DIR              | EXPECTED
        false   | false       | s3('/foo/work/dir')   | 's3://foo/work/dir'
        true    | false       | s3('/foo/work/dir')   | '/fusion/s3/foo/work/dir'
        false   | true        | Path.of('/nfs/work')  | '/nfs/work'
    }

    def 'should get array launch command' (){
        given:
        def executor = Spy(AwsBatchExecutor) {
            isFusionEnabled()>>FUSION
            isWorkDirDefaultFS()>>DEFAULT_FS
            getAwsOptions() >> Mock(AwsOptions) {
                getS5cmdPath() >> { S5CMD ? 's5cmd' : null }
                getAwsCli() >> { 'aws' }
            }
        }
        expect:
        executor.getArrayLaunchCommand(TASK_DIR) == EXPECTED

        where:
        FUSION  | DEFAULT_FS  | S5CMD   | TASK_DIR            | EXPECTED
        false   | false       | false   | 's3://foo/work/dir' | 'bash -o pipefail -c \'trap "{ ret=$?; aws s3 cp --only-show-errors .command.log s3://foo/work/dir/.command.log||true; exit $ret; }" EXIT; aws s3 cp --only-show-errors s3://foo/work/dir/.command.run - | bash 2>&1 | tee .command.log\''
        false   | false       | true    | 's3://foo/work/dir' | 'bash -o pipefail -c \'trap "{ ret=$?; s5cmd cp .command.log s3://foo/work/dir/.command.log||true; exit $ret; }" EXIT; s5cmd cat s3://foo/work/dir/.command.run | bash 2>&1 | tee .command.log\''
        and:
        true    | false       | false   | '/fusion/work/dir'  | 'bash /fusion/work/dir/.command.run'
        false   | true        | false   | '/nfs/work/dir'     | 'bash /nfs/work/dir/.command.run 2>&1 > /nfs/work/dir/.command.log'
    }

}
