/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.trace

import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class CleanupStrategyTest extends Specification {

    def "should check if it's a valid cleanup strategy" () {

        expect:
        CleanupStrategy.isValid(NAME) == EXPECTED

        where:
        NAME            | EXPECTED
        null            | false
        ''              | false
        'foo'           | false
        and:
        'lazy'          | true
        'eager'         | true
        'aggressive'    | true
        and:
        'LAZY'          | true
        'EAGER'         | true
        'AGGRESSIVE'    | true
    }

}
