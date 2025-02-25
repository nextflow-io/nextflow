/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import java.nio.file.Path

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import nextflow.Session
import nextflow.SysEnv
import nextflow.processor.TaskHandler
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleBatchExecutorTest extends Specification {

    def 'should check is fusion' () {
        given:
        SysEnv.push(ENV)
        and:
        def sess = Mock(Session) {
            getConfig() >> CONFIG
        }
        def executor = new GoogleBatchExecutor(session: sess)

        expect:
        executor.isFusionEnabled() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        CONFIG                      | ENV                       | EXPECTED
        [:]                         | [:]                       | false
        [fusion:[enabled: true]]    | [:]                       | true
        [fusion:[enabled: false]]   | [FUSION_ENABLED:'true']   | false     // <-- config has priority
        [:]                         | [FUSION_ENABLED:'true']   | true
        [:]                         | [FUSION_ENABLED:'false']  | false
    }

    @Unroll
    def 'should check cloudinfo enabled' () {
        given:
        SysEnv.push(ENV)
        and:
        def sess = Mock(Session) { getConfig() >> [:] }
        def executor = new GoogleBatchExecutor(session: sess)

        expect:
        executor.isCloudinfoEnabled() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        ENV                             | EXPECTED
        [:]                             | true
        [NXF_CLOUDINFO_ENABLED:'true']  | true
        [NXF_CLOUDINFO_ENABLED:'false'] | false
    }

    def 'should get array index variable and start' () {
        given:
        def executor = Spy(GoogleBatchExecutor)
        expect:
        executor.getArrayIndexName() == 'BATCH_TASK_INDEX'
        executor.getArrayIndexStart() == 0
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(GoogleBatchExecutor)
        expect:
        executor.getArrayTaskId(JOB_ID, TASK_INDEX) == EXPECTED

        where:
        JOB_ID      | TASK_INDEX    | EXPECTED
        'foo'       | 1             | '1'
        'bar'       | 2             | '2'
    }

    protected Path gs(String str) {
        def b = str.tokenize('/')[0]
        def p = str.tokenize('/')[1..-1].join('/')
        CloudStorageFileSystem.forBucket(b).getPath('/'+p)
    }

    @Unroll
    def 'should get array task id' () {
        given:
        def executor = Spy(GoogleBatchExecutor) {
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
        false   | false       | gs('/foo/work/dir')   | '/mnt/disks/foo/work/dir'
        true    | false       | gs('/foo/work/dir')   | '/fusion/gs/foo/work/dir'
        false   | true        | Path.of('/nfs/work')  | '/nfs/work'
    }

    def 'should get array launch command' (){
        given:
        def executor = Spy(GoogleBatchExecutor) {
            isFusionEnabled()>>FUSION
            isWorkDirDefaultFS()>>DEFAULT_FS
        }
        expect:
        executor.getArrayLaunchCommand(TASK_DIR) == EXPECTED

        where:
        FUSION  | DEFAULT_FS  | TASK_DIR            | EXPECTED
        false   | false       | 'gs://foo/work/dir' | '/bin/bash -o pipefail -c \'trap "{ cp .command.log gs://foo/work/dir/.command.log; }" ERR; /bin/bash gs://foo/work/dir/.command.run 2>&1 | tee .command.log\''
        true    | false       | '/fusion/work/dir'  | 'bash /fusion/work/dir/.command.run'
        false   | true        | '/nfs/work/dir'     | 'bash /nfs/work/dir/.command.run 2>&1 > /nfs/work/dir/.command.log'
    }

    def 'should validate shouldDeleteJob method' () {
        given:
        def executor = Spy(GoogleBatchExecutor)

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
}
