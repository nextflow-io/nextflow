/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch

import nextflow.Session
import nextflow.SysEnv
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
}
