/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.executor

import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.aws.batch.AwsBatchExecutor
import nextflow.util.ThrottlingExecutor
import spock.lang.Specification

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

    def 'should kill tasks' () {
        given:
        def reaper = Mock(ThrottlingExecutor) {
            submit(_) >> { args -> args[0]() }
        }
        def executor = Spy(AwsBatchExecutor)
        executor.@reaper = reaper

        when:
        executor.killTask('job-id')
        executor.killTask('job-id')
        then:
        1 * executor.killTask0('job-id') >> null

        when:
        executor.killTask('array-job-id:0')
        executor.killTask('array-job-id:1')
        then:
        1 * executor.killTask0('array-job-id') >> null
    }

}
