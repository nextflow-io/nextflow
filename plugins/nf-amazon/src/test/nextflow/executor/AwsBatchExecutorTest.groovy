/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.executor

import nextflow.Session
import nextflow.cloud.aws.batch.AwsBatchExecutor
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AwsBatchExecutorTest extends Specification {

    def 'should check is fusion' () {
        given:
        def sess = Mock(Session) {
            getConfig() >> CONFIG
        }
        def executor = new AwsBatchExecutor(session: sess, sysEnv: ENV)

        expect:
        executor.isFusionEnabled() == EXPECTED

        where:
        CONFIG                      | ENV                           | EXPECTED
        [:]                         | [:]                           | false
        [fusion:[enabled: true]]    | [:]                           | true
        [fusion:[enabled: false]]   | [NXF_FUSION_ENABLED:'true']   | false     // <-- config has priority
        [:]                         | [NXF_FUSION_ENABLED:'true']   | true
        [:]                         | [NXF_FUSION_ENABLED:'false']  | false

    }

}
