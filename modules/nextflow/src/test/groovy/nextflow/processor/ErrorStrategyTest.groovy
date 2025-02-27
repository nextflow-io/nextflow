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

package nextflow.processor

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ErrorStrategyTest extends Specification {

    def "should check if it's a valid error strategy" () {

        expect:
        ErrorStrategy.isValid(NAME) == EXPECTED

        where:
        NAME            | EXPECTED
        null            | false
        ''              | false
        'foo'           | false
        and:
        'terminate'     | true
        'ignore'        | true
        'retry'         | true
        'finish'        | true
        and:
        'TERMINATE'     | true
        'IGNORE'        | true
        'RETRY'         | true
        'FINISH'        | true
    }

}
