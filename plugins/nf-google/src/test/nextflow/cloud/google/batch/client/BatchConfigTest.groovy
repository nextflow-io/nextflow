/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.cloud.google.batch.client

import nextflow.Session
import nextflow.util.MemoryUnit
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchConfigTest extends Specification {

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should create batch config' () {
        given:
        def CONFIG = [:]
        def session = Mock(Session) { getConfig()>>CONFIG }

        when:
        def config = BatchConfig.create(session)
        then:
        !config.getSpot()
        and:
        config.retryConfig.maxAttempts == 5
        config.maxSpotAttempts == 0
        config.autoRetryExitCodes == [50001]
        and:
        !config.bootDiskImage
        !config.bootDiskSize
    }

    @Requires({System.getenv('GOOGLE_APPLICATION_CREDENTIALS')})
    def 'should create batch config with custom settings' () {
        given:
        def CONFIG = [google: [
            batch: [
                spot: true,
                maxSpotAttempts: 8,
                autoRetryExitCodes: [50001, 50003, 50005],
                retryPolicy: [maxAttempts: 10],
                bootDiskImage: 'batch-foo',
                bootDiskSize: '100GB'
            ]
        ] ]
        def session = Mock(Session) { getConfig()>>CONFIG }

        when:
        def config = BatchConfig.create(session)
        then:
        config.getSpot()
        and:
        config.retryConfig.maxAttempts == 10
        config.maxSpotAttempts == 8
        config.autoRetryExitCodes == [50001, 50003, 50005]
        and:
        config.bootDiskImage == 'batch-foo'
        config.bootDiskSize == MemoryUnit.of('100GB')
    }

}
