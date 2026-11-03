/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.azure.config

import nextflow.Global
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchOptsTest extends Specification {

    @Unroll
    def 'should fetch location' ()  {

        expect:
        new AzBatchOpts(endpoint: ENDPOINT, location: LOC, [:]).getLocation() == EXPECTED

        where:
        LOC                 | ENDPOINT                                  |  EXPECTED
        null                | null                                      | null
        'europe1'           | null                                      | 'europe1'
        'europe1'           | 'https://foo.xyz.batch.azure.com'         | 'europe1'
        null                | 'https://foo.india2.batch.azure.com'      | 'india2'

    }

    @Unroll
    def 'should fetch account name' ()  {

        expect:
        new AzBatchOpts(endpoint: ENDPOINT, accountName: NAME, [:]).getAccountName() == EXPECTED

        where:
        NAME                | ENDPOINT                                  |  EXPECTED
        null                | null                                      | null
        'foo'               | null                                      | 'foo'
        'foo'               | 'https://bar.eu1.batch.azure.com'         | 'foo'
        null                | 'https://bar.eu2.batch.azure.com'         | 'bar'

    }

    def 'should get account name & key'() {
        when:
        def opts1 = new AzBatchOpts([:], [:])
        then:
        opts1.accountName == null
        opts1.accountKey == null

        when:
        def opts2 = new AzBatchOpts(
                [accountName: 'xyz', accountKey: '123'],
                [AZURE_BATCH_ACCOUNT_NAME: 'env-name', AZURE_BATCH_ACCOUNT_KEY:'env-key'])
        then:
        opts2.accountName == 'xyz'
        opts2.accountKey == '123'


        when:
        def opts3 = new AzBatchOpts(
                [:],
                [AZURE_BATCH_ACCOUNT_NAME: 'env-name', AZURE_BATCH_ACCOUNT_KEY:'env-key'])
        then:
        opts3.accountName == 'env-name'
        opts3.accountKey == 'env-key'
    }

    @Unroll
    def 'should check azcopy install' () {
        given:
        Global.session = Mock(Session) { getConfig()>>[fusion:[enabled: FUSION]] }
        AzBatchOpts opts = Spy(AzBatchOpts, constructorArgs: [ CONFIG,[:] ])

        expect:
        opts.getCopyToolInstallMode() == EXPECTED

        where:
        EXPECTED                    | CONFIG                                                    | FUSION
        CopyToolInstallMode.task    | [:]                                                       | false
        CopyToolInstallMode.node    | [allowPoolCreation: true]                                 | false
        CopyToolInstallMode.node    | [autoPoolMode: true]                                      | false
        CopyToolInstallMode.node    | [allowPoolCreation: true, copyToolInstallMode: 'node']    | false
        CopyToolInstallMode.task    | [allowPoolCreation: true, copyToolInstallMode: 'task']    | false
        CopyToolInstallMode.task    | [copyToolInstallMode: 'task']                             | false
        CopyToolInstallMode.node    | [copyToolInstallMode: 'node']                             | false
        CopyToolInstallMode.off     | [copyToolInstallMode: 'off']                              | false
        CopyToolInstallMode.off     | [:]                                                       | true

    }

    def 'should set jobMaxWallClockTime' () {
        when:
        def opts1 = new AzBatchOpts([:], [:])
        then:
        opts1.jobMaxWallClockTime.toString() == '30d'

        when:
        def opts2 = new AzBatchOpts([jobMaxWallClockTime: '3d'], [:])
        then:
        opts2.jobMaxWallClockTime.toString() == '3d'

        when:
        def opts3 = new AzBatchOpts([jobMaxWallClockTime: '12h'], [:])
        then:
        opts3.jobMaxWallClockTime.toString() == '12h'
    }
}
