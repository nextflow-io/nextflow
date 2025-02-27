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

package nextflow.cloud.google.config

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GoogleRetryOptsTest extends Specification {

    def 'should get retry opts' () {
        when:
        def opts1 = new GoogleRetryOpts([:])
        then:
        opts1.maxAttempts == 10
        opts1.multiplier == 2.0d
        opts1.maxDelay ==  Duration.of('90s')

        when:
        def opts2 = new GoogleRetryOpts([maxAttempts: 5, maxDelay: '5s', multiplier: 10])
        then:
        opts2.maxAttempts == 5
        opts2.multiplier == 10d
        opts2.maxDelay ==  Duration.of('5s')
    }
}
