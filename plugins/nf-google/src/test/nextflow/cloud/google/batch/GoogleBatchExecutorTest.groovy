/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.google.batch.client.BatchClient
import spock.lang.Specification
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

    def 'should kill tasks' () {
        given:
        def client = Mock(BatchClient)
        def executor = new GoogleBatchExecutor(client: client)

        when:
        executor.killTask('job-id')
        executor.killTask('job-id')
        then:
        1 * client.deleteJob('job-id')
    }

}
