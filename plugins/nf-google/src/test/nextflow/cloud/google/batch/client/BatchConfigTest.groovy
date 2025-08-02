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
import nextflow.SysEnv
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchConfigTest extends Specification {

    def 'should create batch config' () {
        given:
        SysEnv.push([:])
        and:
        def CONFIG = [google: [
                              project:'proj-123',
                                batch: [
                                    spot: true,
                                    retryPolicy: [maxAttempts: 10]
                                ]
                        ] ]
        def session = Mock(Session) { getConfig()>>CONFIG }

        when:
        def config = BatchConfig.create(session)
        then:
        config.googleOpts.projectId == 'proj-123'
        and:
        config.getSpot()
        and:
        config.retryConfig.maxAttempts == 10
        
        cleanup:
        SysEnv.pop()
    }

}
